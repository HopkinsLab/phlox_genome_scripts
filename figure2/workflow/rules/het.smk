###################
# VCF het #
###################

# Compute FIS values
rule run_vcfhet:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.vcf.gz",
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/vcfhet/populations.snps.lFilt.iFilt.het"
    params:
        wsize = config["pixyParams"]["windowSize"],
        minsize = config["pixyParams"]["minChromSize"]
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 30,
        tmpdir = "/scratch",
        slurm_partition = "serial_requeue"
    threads: 1
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        out_prefix={output}
        out_prefix=${{out_prefix%.*}}
        vcftools --gzvcf {input} --het --temp $TMPDIR --out $out_prefix
        """


