rule bcf_stats:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.vcf.gz"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.vcf.gz.stats"
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * 6000,
        runtime = 720,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    threads: 4
    shell:
        """
        bcftools stats --threads {threads} -s "-" {input} > {output}
        """
