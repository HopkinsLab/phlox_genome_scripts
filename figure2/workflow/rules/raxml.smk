rule prep_raxml:
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
        bcftools view -r chrlg1_RagTag,chrlg2_RagTag,chrlg3_RagTag,chrlg4_RagTag,chrlg5_RagTag,chrlg6_RagTag,chrlg7_RagTag $vcf_fn | \
            bcftools +fill-tags /dev/stdin -- -t AN | \
            bcftools query -i 'AN > 2' -Hf "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n" | \
            workflow/scripts/raxml.prep_seq.pl > $out_phy
        
        """
