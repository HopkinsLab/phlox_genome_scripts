rule mpileup:
    input:
        bam = lambda wildcards: f"{BAM_DIR}/{wildcards.samp}.vs_phlox_pilo.v1.bam"
    output:
        vcf = "results/pileups/{samp}/{samp}.chunk{chunk}.mpileup.vcf.gz",
        idx = "results/pileups/{samp}/{samp}.chunk{chunk}.mpileup.vcf.gz.csi"
    resources:
        mem_mb = 2000,
        runtime = 240,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    threads: 4
    params:
        ref = REF_FASTA,
        chunk_dir = CHUNK_DIR,
        chunk_prefix = CHUNK_PREFIX
    shell:
        """
        export OMP_NUM_THREADS={threads}
        
        chunk_bed={params.chunk_dir}/{params.chunk_prefix}.n{wildcards.chunk}.bed
        tmprun=$(mktemp)
        nthreads=$(expr {threads} - 1)
        echo "samtools view -L $chunk_bed -b {input.bam} | \\" > $tmprun
        echo "bcftools mpileup -X ont -Oz --threads $nthreads -f {params.ref} -a AD,DP /dev/stdin  > $out_fn" >> $tmprun
        srun -c {threads} bash $tmprun

        echo "Indexing vcf"
        tabix -fC {output.vcf}
        """

rule sample_genotypes_from_mpileup:
    input:
        vcf = "results/pileups/{samp}/{samp}.chunk{chunk}.mpileup.vcf.gz",
        idx = "results/pileups/{samp}/{samp}.chunk{chunk}.mpileup.vcf.gz.csi"
    output:
        vcf = "results/sampled_genotypes/{samp}/{samp}.chunk{chunk}.mpileup.rep{repi}.vcf.gz",
        idx = "results/sampled_genotypes/{samp}/{samp}.chunk{chunk}.mpileup.rep{repi}.vcf.gz.csi"
    resources:
        mem_mb = 3000,
        runtime = 720,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    threads: 1
    params:
        chunk_dir = CHUNK_DIR,
        chunk_prefix = CHUNK_PREFIX
    shell:
        """
        export OMP_NUM_THREADS={threads}
        
        chunk_bed={params.chunk_dir}/{params.chunk_prefix}.n{wildcards.chunk}.bed
        echo "Writing output to:"
        echo "$out_file"
        bcftools view -R $chunk_bed -i 'FORMAT/DP > 0' {input.vcf} | \
            workflow/scripts/sample_genotypes_from_mpileup.pl --out_file {output.vcf}
        tabix -fC {output.vcf}
        """

rule filter_sampled_gts:
    input:
        vcf = "results/sampled_genotypes/{samp}/{samp}.chunk{chunk}.mpileup.rep{repi}.vcf.gz",
        idx = "results/sampled_genotypes/{samp}/{samp}.chunk{chunk}.mpileup.rep{repi}.vcf.gz.csi"
    output:
        vcf = "results/sampled_genotypes/filtered/chunk{chunk}/{samp}.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.vcf.gz",
        idx = "results/sampled_genotypes/filtered/chunk{chunk}/{samp}.chunk{chunk}.mpileup.rep{repi}.DP_{min_dp}.vcf.gz.csi"
    resources:
        mem_mb = 2000,
        runtime = 120,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    threads: 3
    shell:
        """
        srun -c {threads} bcftools view -o {output.vcf} --threads {threads} -Oz -i "FORMAT/DP >= {wildcards.min_dp}" {intput.vcf}
        tabix -fC {output.vcf}
        """
