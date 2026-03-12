################################################################
# Stacks

rule prep_gstacks:
    input: 
        config["metadata"]
    output: 
        popmap = "popmaps/{breadth}/popmap.{sp_grp}.txt",
        donefile = touch("results/{breadth}/{sp_grp}/res_align.done")
    params:
        bam_dir = config["bamDir"],
        species = lambda wildcards: config["species"][wildcards.breadth][wildcards.sp_grp],
        rm_samps = config["removeSamples"],
        awkscript = "workflow/scripts/make_popmap.awk"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["prep_gstacks"]["mem_mb"],
        runtime = res_config["prep_gstacks"]["runtime"],
        tmpdir = res_config["prep_gstacks"]["tmpdir"],
        partition = res_config["prep_gstacks"]["partition"]
    shell:
        """
        # Make popmap file
        sp_arr=( {params.species} )
        if [[ ${{#sp_arr[@]}} -eq 1 ]]; then
            if [[ ${{sp_arr[0]}} == "pilo" ]]; then
                pop_col="species"
            else
                pop_col="Population"
            fi
        else
            pop_col="species"
        fi

        for sp in ${{sp_arr[@]}} ; do
            awk -v sp=$sp -v rm_samps={params.rm_samps} -v pop_col=$pop_col -f {params.awkscript} {input} >> {output.popmap}
        done

        # Link bams
        bam_dir={output.donefile}
        bam_dir=${{bam_dir%.done}}
        mkdir -p $bam_dir

        for samp in $(cut -f 1 {output.popmap}); do
            cp {params.bam_dir}/$samp.bam {params.bam_dir}/$samp.bam.csi $bam_dir
        done
        """


# Stacks pipeline: Run gstacks
rule run_gstacks:
    input: 
        popmap = "popmaps/{breadth}/popmap.{sp_grp}.txt",
        donefile = "results/{breadth}/{sp_grp}/res_align.done"
    output:
        "results/{breadth}/{sp_grp}/res_stacks/catalog.fa.gz"
    params:
        stacks_dir = config["stacksDir"]
    threads: res_config["run_gstacks"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["run_gstacks"]["mem_mb"],
        runtime = res_config["run_gstacks"]["runtime"],
        tmpdir = res_config["run_gstacks"]["tmpdir"],
        partition = res_config["run_gstacks"]["partition"]
    shadow: 
        config["shadowLevel"]
    shell:
        """
        module load gcc/12.2.0-fasrc01
        
        bam_dir={input.donefile}
        bam_dir=${{bam_dir%.done}}

        {params.stacks_dir}/gstacks -I $bam_dir -M {input.popmap} -O $( dirname {output} ) -t {threads} --rm-unpaired-reads
        """

# Stacks pipeline: Run populations
rule run_populations:
    input:
        popmap = "popmaps/{breadth}/popmap.{sp_grp}.txt",
        catalog = "results/{breadth}/{sp_grp}/res_stacks/catalog.fa.gz"
    output:
        snps_vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}/populations.snps.vcf",
        all_vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}/populations.all.vcf"
    params:
        stacks_dir = config["stacksDir"],
    threads: res_config["run_populations"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["run_populations"]["mem_mb"],
        runtime = res_config["run_populations"]["runtime"],
        tmpdir = res_config["run_populations"]["tmpdir"],
        partition = res_config["run_populations"]["partition"]
    shadow: 
        config["shadowLevel"]
    shell:
        """
        module load gcc/12.2.0-fasrc01
        {params.stacks_dir}/populations \
            -M {input.popmap} \
            -P $( dirname {input.catalog} ) \
            -O $( dirname {output.snps_vcf} ) \
            --min-samples-per-pop {wildcards.min_samples_per_pop} \
            --min-samples-overall {wildcards.min_samples_overall} \
            --vcf --vcf-all -t {threads}
        """
################################################################
# Filtering and Plink

rule filter_vcf:
    input:
        vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/populations.{vcf_type}.vcf",
        popmap = "popmaps/{breadth}/popmap.{sp_grp}.txt",
        metadata = "metadata.tsv"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/populations.{vcf_type}.lFilt.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["filter_vcf"]["mem_mb"],
        runtime = res_config["filter_vcf"]["runtime"],
        tmpdir = res_config["filter_vcf"]["tmpdir"],
        partition = res_config["filter_vcf"]["partition"]
    # shadow:
    #     config["shadowLevel"]
    threads:
        res_config["filter_vcf"]["threads"]
    params:
        plink_dir = config["plinkDir"],
        min_DP = config["filterSettings"]["min_DP"],
        max_DP = config["filterSettings"]["max_DP"],
        min_GQ = config["filterSettings"]["min_GQ"]
    shell:
        """
        if [[ ! -d $TMPDIR ]]; then
            mkdir $TMPDIR
        fi

        awk -v popmap={input.popmap} \
            -v min_samp_per_pop={wildcards.min_samples_per_pop} \
            -v min_samp_overall={wildcards.min_samples_overall} \
            -v min_DP={params.min_DP} \
            -v max_DP={params.max_DP} \
            -v min_GQ={params.min_GQ} \
            -f workflow/scripts/filter_variants.awk \
            {input.vcf} | \
        bcftools sort -Oz > {output}
        
        tabix -C {output}
        
        # Link csi
        cd $(dirname {output})
        ln -fs $(basename {output}).csi $(basename {output}).tbi
        """

rule vcf_to_plink:
    input:
        snps_vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/populations.snps.lFilt.vcf.gz",
        all_vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/populations.all.lFilt.vcf.gz",
        popmap = "popmaps/{breadth}/popmap.{sp_grp}.txt",
    output:
        bim = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.bim",
        snps_vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.vcf.gz",
        all_vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.vcf.gz"
    params:
        plink_dir = config["plinkDir"]
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["vcf_to_plink"]["mem_mb"],
        runtime = res_config["vcf_to_plink"]["runtime"],
        tmpdir = res_config["vcf_to_plink"]["tmpdir"],
        partition = res_config["vcf_to_plink"]["partition"]
    shell:
        """
        out_prefix={output.bim}
        out_prefix=${{out_prefix%.*}}
        
        tmp_dir=$(mktemp -d)

        tmp_prefix=$tmp_dir/plink
        tmp_samp=$tmp_prefix.samps.txt

        # Make sample update file
        module load R/4.2.2-fasrc01
        workflow/scripts/plink.update_samples.R metadata.tsv {input.popmap} > $tmp_samp


        bcftools view {input.snps_vcf} | \
            {params.plink_dir}/plink --vcf /dev/stdin \
                --double-id --allow-extra-chr \
                --make-bed --out $tmp_prefix

        {params.plink_dir}/plink --bfile $tmp_prefix \
            --allow-extra-chr \
            --update-ids $tmp_samp \
            --mind {wildcards.mind} \
            --make-bed --out $out_prefix
        

        ########################################################################
        tmp_samp=$tmp_prefix.samps.keep.txt
        cut -f 2 -d ' ' $out_prefix.fam > $tmp_samp

        bcftools view -Oz --no-version -S $tmp_samp {input.snps_vcf} > {output.snps_vcf}
        tabix -fC {output.snps_vcf}
        
        bcftools view -Oz --no-version -S $tmp_samp {input.all_vcf}  > {output.all_vcf}
        tabix -fC {output.all_vcf}

        cd $(dirname {output.snps_vcf})
        rm -f $(basename {output.snps_vcf}).tbi
        ln -fs $(basename {output.snps_vcf}).csi $(basename {output.snps_vcf}).tbi
        rm -f $(basename {output.all_vcf}).tbi
        ln -fs $(basename {output.all_vcf}).csi  $(basename {output.all_vcf}).tbi
        
        #### Cleanup
        rm -rf $tmp_prefix*
        """
rule make_popfiles:
    input:
        vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.vcf.gz",
        popmap = "popmaps/{breadth}/popmap.{sp_grp}.txt"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.pop"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["make_popfiles"]["mem_mb"],
        runtime = res_config["make_popfiles"]["runtime"],
        tmpdir = res_config["make_popfiles"]["tmpdir"],
        partition = res_config["make_popfiles"]["partition"]
    shell:
        """
        # make populations file
        > {output}
        while read samp; do
            awk -v samp=$samp '$1 == samp {{print $0; exit 0}}' {input.popmap} >> {output}
        done < <( zcat {input.vcf} | grep -m 1 '^#CHROM' | cut -f 10- | tr '\t' '\n' )
        """
