##############
# PCA IMPUTE #
##############
rule prep_impute:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.bim"
    output:
        traw = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca_impute/populations.snps.lFilt.iFilt.pca.impute.traw"
    params:
        plink_dir = config["plinkDir"],
        ld_prune = config["pcaImputeParams"]["ldPrune"],
        hwe      = config["pcaImputeParams"]["hwe"],
        maf      = config["pcaImputeParams"]["maf"]
    resources:
        mem_mb_per_cpu = 4000,
        runtime = 20,
        tmpdir = "/scratch",
        slurm_partition = "serial_requeue"
    threads: 1
    shell:
        """
        in_prefix={input}
        in_prefix=${{in_prefix%.*}}

        out_prefix={output.traw}
        out_prefix=${{out_prefix%.impute*}}

        # LD Prune
        {params.plink_dir}/plink --bfile $in_prefix \
            --allow-extra-chr \
            --indep-pairwise {params.ld_prune} \
            --out $out_prefix
        mv -f $out_prefix.log $out_prefix.prune.log

        # Export transposed because we need the tfam file
        {params.plink_dir}/plink --bfile $in_prefix \
            --allow-extra-chr \
            --extract $out_prefix.prune.in \
            --hwe {params.hwe} \
            --maf 0.025 \
            --recode 12 transpose \
            --out $out_prefix.impute
        mv -f $out_prefix.impute.log $out_prefix.impute.transpose.log

        # Export in R-friendly traw format
        {params.plink_dir}/plink --tfile $out_prefix.impute \
            --allow-extra-chr \
            --recode 12 A-transpose \
            --out $out_prefix.impute
        mv -f $out_prefix.impute.log $out_prefix.impute.A-transpose.log
        """

rule run_impute:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca_impute/populations.snps.lFilt.iFilt.pca.impute.traw"
    output:
        tped   = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca_impute/populations.snps.lFilt.iFilt.pca.impute.missForest.tped",
        err_fn = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca_impute/populations.snps.lFilt.iFilt.pca.impute.missForest.OOBerror.txt"
    params:
        maxiter = config["pcaImputeParams"]["maxiter"]
    resources:
        mem_mb_per_cpu = 2000,
        runtime = 720,
        tmpdir = "/scratch",
        slurm_partition = "serial_requeue"
    threads: 8
    shell:
        """
        # Impute
        echo ""
        echo "$(date): Imputing data"
        module load R/4.2.2-fasrc01

        in_prefix={input}
        in_prefix=${{in_prefix%.*}}
        
        out_prefix={output.tped}
        out_prefix=${{out_prefix%.*}}

        use_seed=$RANDOM
        echo $use_seed > $out_prefix.seed
        
        workflow/scripts/impute.R $in_prefix $out_prefix {params.maxiter} $SLURM_CPUS_PER_TASK $use_seed
        """


rule pca_impute:
    input:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca_impute/populations.snps.lFilt.iFilt.pca.impute.missForest.tped",
    output:
        eigval = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca_impute/populations.snps.lFilt.iFilt.pca.impute.missForest.eigenval",
        eigvec = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca_impute/populations.snps.lFilt.iFilt.pca.impute.missForest.eigenvec"
    params:
        plink_dir = config["plinkDir"],
    resources:
        mem_mb_per_cpu = 3000,
        runtime = 30,
        tmpdir = "/scratch",
        slurm_partition = "serial_requeue"
    threads: 1
    shell:
        """
        in_prefix={input}
        in_prefix=${{in_prefix%.*}}

        out_prefix={output.eigvec}
        out_prefix=${{out_prefix%.*}}

        # PCA
        echo ""
        echo "PCA..."
        {params.plink_dir}/plink \
            --tfile $in_prefix \
            --allow-extra-chr \
            --pca \
            --out $out_prefix
        """
