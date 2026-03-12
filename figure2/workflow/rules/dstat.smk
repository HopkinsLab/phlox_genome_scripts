rule dstat_quartets:
    input:
        vcf = "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.vcf.gz",
        ref = config["ref"],
        popmap = "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.pop",
    output:
        "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/dstat/populations.snps.lFilt.iFilt.patt_counts.tsv"
    resources:
        mem_mb = 16000,
        runtime = 120,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    params:
        nreps = 3 # How many replicates each drummondii individual should get
    shell:
        """
        module load R/4.2.2-fasrc01

        bcftools view -M 2 -i 'TYPE="snp"' {input.vcf} | \
            awk -f workflow/scripts/flip_strand.ddrad.awk | \
            bcftools norm -c s -f {input.ref} 2> {output}.norm.log | \
            bcftools query -Hf '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | \
            sed -e '1 s/^#//' -e '1 s/\[[0-9]\+\]//g'  -e '1 s/:GT//g' -e 's/0\/0/0/g' -e 's/1\/1/1/g' | \
            workflow/scripts/dstat.quartets.R {input.popmap} {params.nreps} {output}
        """

rule plot_dstat_quartets:
    input:
        tsv = "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/dstat/populations.snps.lFilt.iFilt.patt_counts.tsv",
        metadata = config["metadata"]
    output:
        touch("results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/dstat/populations.snps.lFilt.iFilt.dstat.plots.done")
    resources:
        mem_mb = 2000,
        runtime = 60,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        module load R/4.2.2-fasrc01

        out_prefix={output}
        out_prefix=${{out_prefix%.*}}

        workflow/scripts/plot_dstat.R {input.tsv} {input.metadata} $out_prefix
        """
