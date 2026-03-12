#######
# PCA #
#######

# Run a plink PCA with some pre-filtering
rule run_pca:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.bim"
    output:
        eigenval = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca/populations.snps.lFilt.iFilt.pca.eigenval",
        eigenvec = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca/populations.snps.lFilt.iFilt.pca.eigenvec"
    params:
        plink_dir = config["plinkDir"],
        ld_prune = config["pcaParams"]["ldPrune"],
        hwe      = config["pcaParams"]["hwe"]
    resources:
        mem_mb = 4000,
        runtime = 60,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        in_prefix={input}
        in_prefix=${{in_prefix%.*}}

        out_prefix={output.eigenval}
        out_prefix=${{out_prefix%.*}}


        # LD Prune
        {params.plink_dir}/plink --bfile $in_prefix \
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

        {params.plink_dir}/plink --bfile $out_prefix \
            --allow-extra-chr \
            --pca --out $out_prefix
        """

rule plot_pca:
    input:
        eigenval = "results/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca/populations.snps.lFilt.iFilt.pca.eigenval",
        eigenvec = "results/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca/populations.snps.lFilt.iFilt.pca.eigenvec",
        metadata = "metadata.tsv"
    output:
        touch("results/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca/plots/populations.snps.lFilt.iFilt.pca.plots.done")
    resources:
        mem_mb = 4000,
        runtime = 45,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    shell:
        """
        in_prefix={input.eigenval}
        in_prefix=${{in_prefix%.*}}

        out_prefix={output}
        out_prefix=${{out_prefix%.*}}

        module load gcc/12.2.0-fasrc01
        module load R/4.2.2-fasrc01
        workflow/scripts/plot_pca.R $in_prefix {input.metadata} $out_prefix
        """
