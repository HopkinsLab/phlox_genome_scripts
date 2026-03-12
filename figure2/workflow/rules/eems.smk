########
# EEMS #
########

# Prepare data for an EEMS run
rule prep_eems:
    input:
        metadata = config["metadata"],
        bim      = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.bim"
    output:
        runfile = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/eems/{sp}/populations.snps.lFilt.iFilt.eems.param",
        diffs   = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/eems/{sp}/populations.snps.lFilt.iFilt.eems.diffs",
        coord   = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/eems/{sp}/populations.snps.lFilt.iFilt.eems.coord",
        outer   = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/eems/{sp}/populations.snps.lFilt.iFilt.eems.outer"
    params:
        plink_dir   = config["plinkDir"],
        bed2diffs   = config["eemsBed2Diffs"],
        habitat_dir = config["eemsHabitatDir"],
        rm_samps    = config["eemsParams"]["rmSampsPath"],
        nDemes      = config["eemsParams"]["nDemes"],
        diploid     = config["eemsParams"]["diploid"],
        numMCMCIter = config["eemsParams"]["numMCMCIter"],
        numBurnIter = config["eemsParams"]["numBurnIter"],
        numThinIter = config["eemsParams"]["numThinIter"],
        ld_prune    = config["eemsParams"]["ldPrune"],
        hwe         = config["eemsParams"]["hwe"]
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 60,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        module load gcc/12.2.0-fasrc01
        OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

        in_prefix={input.bim}
        in_prefix=${{in_prefix%.*}}

        out_prefix={output.runfile}
        out_prefix=${{out_prefix%.*}}

        rm_flag=""
        if [[ ! -z {params.rm_samps} ]]; then
            rm_flag="--remove {params.rm_samps}"
        fi

        tmp_prefix=$TMPDIR/$(echo {input.bim} | md5sum | cut -f 1 -d ' ').$RANDOM
   
        ### Keep only chromosomal lgs, samples in current species group ###
        grep '^chrlg' $in_prefix.bim | cut -f 2 > $out_prefix.chr.in
        > $tmp_prefix.keep
        while read -r fid iid ; do
            keepsamp=$( awk -v sp={wildcards.sp} -v iid="$iid" '$2 == iid {{print $7 == sp}}' {input.metadata} )
            if [ $keepsamp -ne 0 ]; then
                echo "$fid $iid" >> $tmp_prefix.keep
            fi 
        done < <( cut -f 1-2 -d ' ' $in_prefix.fam )

        {params.plink_dir}/plink --bfile $in_prefix \
            --allow-extra-chr \
            --extract $out_prefix.chr.in \
            --keep $tmp_prefix.keep \
            --make-bed --out $tmp_prefix


        ### LD Prune and filter ###
        {params.plink_dir}/plink --bfile $tmp_prefix \
            --allow-extra-chr \
            --indep-pairwise {params.ld_prune} \
            --out $out_prefix
        mv $out_prefix.log $out_prefix.prune.log

        {params.plink_dir}/plink --bfile $tmp_prefix \
            --allow-extra-chr \
            --extract $out_prefix.prune.in \
            $rm_flag \
            --hwe {params.hwe} \
            --make-bed --out $out_prefix
        mv $out_prefix.log $out_prefix.filter.log

        # Convert chromosomes to numbers
        for i in $(seq 7); do
            sed -i 's/^chrlg'$i'_RagTag/'$i'/' $out_prefix.bim
        done


        ### Make .diffs file ###
        {params.bed2diffs} --bfile $out_prefix --nthreads {threads}

        ### Make .coord file ###
        workflow/scripts/coord_from_fam.sh {input.metadata} $out_prefix.fam {output.coord}

        ### Copy .outer file ###
        cp {params.habitat_dir}/annuals.outer {output.outer}

        ### Make runfile ###
        n_indiv=$( wc -l $out_prefix.fam | awk '{{print $1}}' )
        n_sites=$( wc -l $out_prefix.bim | awk '{{print $1}}' )

        echo "datapath = $out_prefix"              > {output.runfile}
        echo "mcmcpath = ${{out_prefix}}_out"     >> {output.runfile}
        echo "nIndiv = $n_indiv"                  >> {output.runfile}
        echo "nSites = $n_sites"                  >> {output.runfile}
        echo "nDemes = {params.nDemes}"           >> {output.runfile}
        echo "diploid = {params.diploid}"         >> {output.runfile}
        echo "numMCMCIter = {params.numMCMCIter}" >> {output.runfile}
        echo "numBurnIter = {params.numBurnIter}" >> {output.runfile}
        echo "numThinIter = {params.numThinIter}" >> {output.runfile}
        """

# Run EEMS
rule run_eems:
    input:        
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/eems/{sp}/populations.snps.lFilt.iFilt.eems.param",
    output:
        touch("results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/eems/{sp}/populations.snps.lFilt.iFilt.eems.done"),
    params:
        runeems = config["eemsRunEems"]
    resources:
        mem_mb = 4000,
        runtime = 720,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        n=$( workflow/scripts/seed_from_str.sh {input} )
        
        {params.runeems} --params {input} --seed $n
        """

#rule plot_eems:
#    input:        
#        "res_eems/{vcfname}.{sampgroup}.chrOnly_{chr_only}.minGQ_{min_GQ}.eems.done"
#    output:
#        touch("res_eems/{vcfname}.{sampgroup}.chrOnly_{chr_only}.minGQ_{min_GQ}.eems_plot/plots.done")
#    resources:
#        mem_mb = 4000,
#        runtime = 60,
#        tmpdir = "/scratch",
#        partition = "serial_requeue"
#    shell:
#        """
#        module load R/4.1.0-fasrc01
#        module load gcc/7.1.0-fasrc01
#        module load geos/3.9.1-fasrc01
#        module load gdal/3.2.2-fasrc01
#
#        in_dir={input}
#        in_dir=${{in_dir%.done}}_out
#
#        out_prefix={output}
#        out_prefix=${{out_prefix%.done}}
#
#        workflow/scripts/plot_eems.R $in_dir $out_prefix
#        tar -cvzf $(dirname $out_prefix).tar.gz $(dirname $out_prefix)
#        """
