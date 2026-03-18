
rule sample_snp_patterns_windowed:
    input:
        vcf = expand(
            "{vcf_dir}/{samp}.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.vcf.gz",
            vcf_dir=VCF_DIR,
            samp=SAMPLES,
            allow_missing=True
        ),
        idx = expand(
            "{vcf_dir}/{samp}.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.vcf.gz.csi",
            vcf_dir=VCF_DIR,
            samp=SAMPLES,
            allow_missing=True
        ),
        popfile = config["popFile"]
    output:
        "results/snp_patt/merged.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.patt_counts.w_250kb.tsv"
    resources:
        mem_mb = 2000 * 9,
        runtime = 60 * 8,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    threads: 9
    shell:
        """
        export OMP_NUM_THREADS={threads}
        
        # Fix a naming inconsistency
        tmppop=$(mktemp)
        sed -e 's/^ppilo/ppil/' {input.popfile} > $tmppop

        tmpsh=$(mktemp)
        > $tmpsh
        echo "bcftools merge -m all {input.vcf} | \\" >> $tmpsh
        echo "bcftools view -i 'AN >= 8' | \\" >> $tmpsh
        echo "    ./sample_snp_patterns.windowed.py \\" >> $tmpsh
        echo "        --populations $tmppop \\" >> $tmpsh
        echo "        --outfile {output} \\" >> $tmpsh
        echo "        --pop_order roem,drum,cusp,pilo \\" >> $tmpsh
        echo "        --n_cores {threads} \\" >> $tmpsh
        echo "        --chunk_length 120000 \\" >> $tmpsh
        echo "        --window_size 250000" >> $tmpsh

        srun -c {threads} bash $tmpsh
        """

rule sample_sv_patterns_windowed:
    input:
        vcf = SV_VCF,
        idx = SV_VCF + ".csi",
        popfile = config["popFile"]
    output:
        "results/svs_patt/svs.rep{repi}.SVLEN_{min_svlen}.patt_counts.w_250kb.tsv"
    resources:
        mem_mb = 2000 * 5,
        runtime = 60 * 8,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    threads: 5
    shell:
        """
        export OMP_NUM_THREADS={threads}
        
        tmpvcf=$(mktemp --suffix=".vcf")
        outgrp=$(grep pilo {input.popfile})
        bcftools view -s ^pdrum_REF -i '(INFO/SVLEN > '{wildcards.min_svlen}' || INFO/SVLEN < -'$min_svlen')' {input.vcf} chrlg1_RagTag chrlg2_RagTag chrlg3_RagTag chrlg4_RagTag chrlg5_RagTag chrlg6_RagTag chrlg7_RagTag | \
            bcftools annotate -x '^FORMAT/GT' | \
            awk -v outgrp=$outgrp -v OFS="\t" '/^#CHROM/ {{$(NF+1) = outgrp}} !/^#/ {{$(NF+1) = "0/0"}} {{print $0}}' \
            > $tmpvcf

        tmpsh=$(mktemp)
        > $tmpsh
        
        echo "./sample_snp_patterns.windowed.py \\" >> $tmpsh
        echo "    --vcf $tmpvcf \\" >> $tmpsh
        echo "    --populations {input.popfile} \\" >> $tmpsh
        echo "    --outfile {output} \\" >> $tmpsh
        echo "    --pop_order roem,drum,cusp,pilo \\" >> $tmpsh
        echo "    --n_cores {threads} \\" >> $tmpsh
        echo "    --chunk_length 120000 \\" >> $tmpsh
        echo "    --window_size 250000" >> $tmpsh

        srun -c {threads} bash $tmpsh
        """



rule missing_data_in_vcfs:
    input:
        vcf = expand(
            "{vcf_dir}/{samp}.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.vcf.gz",
            vcf_dir=VCF_DIR,
            samp=SAMPLES,
            allow_missing=True
        ),
        idx = expand(
            "{vcf_dir}/{samp}.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.vcf.gz.csi",
            vcf_dir=VCF_DIR,
            samp=SAMPLES,
            allow_missing=True
        ),
        popfile = config["popFile"],
    output:
        "results/missingness/merged.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.missingness.tsv"
    resources:
        mem_mb = 2000 * 3,
        runtime = 60 * 4,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    threads: 3
    params:
        chunk_dir = CHUNK_DIR,
        chunk_prefix = CHUNK_PREFIX
    shell:
        """
        export OMP_NUM_THREADS={threads}
        chunk_bed={params.chunk_dir}/{params.chunk_prefix}.n{wildcards.chunk}.bed

        # Fix a naming inconsistency
        tmppop=$(mktemp)
        sed -e 's/^ppilo/ppil/' {input.popfile} > $tmppop

        tmprun=$(mktemp)
        > $tmprun
        echo "bcftools merge --threads {threads} -Ov -R $bed_fn -m all {input.vcf} | \\" >> $tmprun
        echo "    bcftools query -f '%CHROM\t%POS[\t%GT]\n' | \\" >> $tmprun
        echo "    awk -v wsize=250000 -f workflow/scripts/missing_data_in_vcfs.awk > $out_fn" >> $tmprun

        srun -c {threads} bash $tmprun
        """

