#############
# ADMIXTURE #
#############
rule prep_admixture:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.bim"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/admix/populations.snps.lFilt.iFilt.admix.bed"
    params:
        plink_dir = config["plinkDir"],
        ld_prune = config["admixParams"]["ldPrune"],
        # maf      = config["admixParams"]["maf"],
        # mac      = config["admixParams"]["mac"],
        # mind     = config["admixParams"]["mind"],
        # geno     = config["admixParams"]["geno"],
        hwe      = config["admixParams"]["hwe"]
    resources:
        mem_mb = 2000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        in_prefix={input}
        in_prefix=${{in_prefix%.*}}

        out_prefix={output}
        out_prefix=${{out_prefix%.bed}}

        tmp_prefix=$TMPDIR/$(echo {output} | md5sum | cut -f 2 -d ' ')

        # Keep only chromosomal lgs
        grep '^chrlg' $in_prefix.bim | cut -f 2 > $out_prefix.chr.in
        {params.plink_dir}/plink --bfile $in_prefix \
            --allow-extra-chr \
            --extract $out_prefix.chr.in \
            --make-bed --out $tmp_prefix

        # LD Prune
        {params.plink_dir}/plink --bfile $tmp_prefix \
            --allow-extra-chr \
            --indep-pairwise {params.ld_prune} \
            --out $out_prefix
        mv $out_prefix.log $out_prefix.prune.log

        {params.plink_dir}/plink --bfile $in_prefix \
            --allow-extra-chr \
            --extract $out_prefix.prune.in \
            --hwe {params.hwe} \
            --make-bed --out $out_prefix
        mv $out_prefix.log $out_prefix.filter.log

        # Convert chromosomes to numbers
        for i in $(seq 7); do
            sed -i 's/^chrlg'$i'_RagTag/'$i'/' $out_prefix.bim
        done
        
                
        """

rule run_admixture:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/admix/populations.snps.lFilt.iFilt.admix.bed"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/admix/populations.snps.lFilt.iFilt.admix.{kval}.log"
    params:
        admix_dir = config["admixDir"],
        n_reps = config["admixParams"]["nReps"]
    resources:
        mem_mb = 4000,
        runtime = 960,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    threads: 4
    shell:
        """
        out_prefix={output}
        out_prefix=${{out_prefix%.log}}

        wd=$( pwd )

        tmpdir=$TMPDIR/$(echo {output}.$RANDOM | md5sum | cut -f 1 -d ' ')
        mkdir $tmpdir; cd $tmpdir
        tmp_log=$(basename {output})

        max_loglik=""
        for i in $(seq {params.n_reps}); do
            seed=$(expr $RANDOM + $(date +%S))

            {params.admix_dir}/admixture -j{threads} --cv --seed=$seed  \
                $wd/{input} {wildcards.kval} | \
                tee $tmp_log
            
            # Get the log likelihood of this run
            cur_loglik=$(grep '^Loglikelihood' $tmp_log | awk '{{print $2}}')

            
            echo "Rep $i: Seed = $seed ; Current_Loglikelihood = $cur_loglik ; PrevMax_Loglikelihood = $max_loglik"
            if [[ -z $max_loglik ]] || [[ $(echo "$cur_loglik > $max_loglik" | bc -l) -eq 1 ]]; then
                mv -f $tmp_log $wd/{output}
                mv -f *.P $wd/$out_prefix.P
                mv -f *.Q $wd/$out_prefix.Q
                max_loglik=$cur_loglik
            else
                rm -f *.log *.P *.Q
            fi
        done
        """


rule plot_admixture:
    input:
        logfiles = ["results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/admix/" + x for x in expand("populations.snps.lFilt.iFilt.admix.{kval}.log", kval=config["admixParams"]["kVals"])],
        metadata = "metadata.tsv"
    output:
        cv_report = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/admix/populations.snps.lFilt.iFilt.admix.CV_report.txt",
        donefile = touch("results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/admix/plots/populations.snps.lFilt.iFilt.admix.CV_report.plots.done")
    params:
        k_vals = config["admixParams"]["kVals"]
    resources:
        mem_mb = 4000,
        runtime = 45,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        grep -h CV {input.logfiles} > {output.cv_report}
        
        in_prefix={output.cv_report}
        in_prefix=${{in_prefix%.CV_report*}}

        out_prefix={output.donefile}
        out_prefix=${{out_prefix%.*}}

        module load gcc/12.2.0-fasrc01
        module load R/4.2.2-fasrc01
        workflow/scripts/plot_admixture.R $in_prefix $(echo {params.k_vals} | tr ' ' ',') {input.metadata} $out_prefix
        """
