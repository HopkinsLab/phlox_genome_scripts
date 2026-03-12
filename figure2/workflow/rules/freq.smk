rule compile_freq:
    input:
        vcf = "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.vcf.gz",
        popmap = "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.pop",
        ref = config["ref"]
    output:
        "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/freq/populations.all.lFilt.iFilt.freq.tsv"
    resources:
        mem_mb = 1000,
        runtime = 60,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -M 2 -i 'TYPE="snp"' {input.vcf} | \
            awk -f workflow/scripts/flip_strand.ddrad.awk | \
            bcftools norm -c s -f {input.ref} | \
            bcftools +fill-tags /dev/stdin -- -t AN,AC,AC_Hom,AC_Het -S {input.popmap} | \
            bcftools query -Hf '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AN_cusp\t%AN_drum\t%AN_roem\t%AN_pilo\t%AC\t%AC_cusp\t%AC_drum\t%AC_roem\t%AC_pilo\t%AC_Hom\t%AC_Hom_cusp\t%AC_Hom_drum\t%AC_Hom_roem\t%AC_Hom_pilo\t%AC_Het\t%AC_Het_cusp\t%AC_Het_drum\t%AC_Het_roem\t%AC_Het_pilo\n' | \
            sed -e '1 s/^#//' -e '1 s/\[[0-9]\+\]//g' > {output}
        """

rule plot_freq:
    input:
        "results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/freq/populations.all.lFilt.iFilt.freq.tsv"
    output:
        touch("results/multi/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/freq/populations.all.lFilt.iFilt.freq.plots.done")
    resources:
        mem_mb = 4000,
        runtime = 15,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        module load gcc/12.2.0-fasrc01
        module load R/4.2.2-fasrc01

        out_prefix={output}
        out_prefix=${{out_prefix%.*}}

        ./workflow/scripts/plot_freq.R {input} $out_prefix
        """

