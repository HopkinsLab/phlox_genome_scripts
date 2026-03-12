
rule prep_treemix:
    input:
        # vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.snps.lFilt.iFilt.vcf.gz",
        pca_out = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pca/populations.snps.lFilt.iFilt.pca.eigenval",
        popfile = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.pop"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/treemix/populations.snps.lFilt.iFilt.input.gz"
    params:
        plink_dir = config["plinkDir"]
    resources:
        mem_mb = 4000,
        runtime = 60,
        tmpdir = "results/tmp",
        partition = "serial_requeue,shared"
    shell:
        """
        module load gcc/12.2.0-fasrc01 
        module load R/4.2.2-fasrc01
        
        in_prefix={input.pca_out}
        in_prefix=${{in_prefix%.eig*}}

        # Check that sample order is the same in fam and popfile
        hash1=$( cut -f 2 -d ' ' $in_prefix.fam | md5sum | cut -f 1 -d ' ' )
        hash2=$( cut -f 1 {input.popfile} | md5sum | cut -f 1 -d ' ' )
        if [ $hash1 != $hash2 ]; then
            echo "input fam and popfile do not have same sample ordering or number" > /dev/stderr
            echo "fam: $in_prefix.fam" > /dev/stderr
            echo "popfile: {input.popfile}" > /dev/stderr
            exit 1
        fi

        # Temp directory
        tmp_prefix=$TMPDIR/$(echo {input.pca_out} | md5sum | cut -f 1 -d ' ')/prep_treemix
        mkdir -p $(dirname $tmp_prefix)

        paste -d ' ' <( cut -f 1-2 -d ' ' $in_prefix.fam ) <( cut -f 2 {input.popfile} ) > $tmp_prefix.clust

        # Calculate frequencies
        {params.plink_dir}/plink --bfile $in_prefix \
            --allow-extra-chr \
            --within $tmp_prefix.clust \
            --freq \
            --out $tmp_prefix
        
        # Convert output to treemix format
        workflow/scripts/prep_treemix.R $tmp_prefix.frq.strat | gzip > {output}
        """
        
        
rule run_treemix:
    input:
        snps = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/treemix/populations.snps.lFilt.iFilt.input.gz",
        # tree = "trees/{breadth}/{sp_grp}.nwk"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/treemix/populations.snps.lFilt.iFilt.input.cov.gz"
    params:
        treemix_dir = config["treemixDir"],
        m = config["treemixParams"]["m"],
    resources:
        mem_mb = 4000,
        runtime = 120,
        tmpdir = "results/tmp",
        partition = "serial_requeue,shared"
    shell:
        """
        module load gcc/12.2.0-fasrc01 

        out_prefix={output}
        out_prefix=${{out_prefix%.cov*}}
        
        {params.treemix_dir}/treemix -i {input.snps} \
            -m {params.m} \
            -root pilo \
            -o $out_prefix
        """

# rule plot_treemix:
#     input:
#         "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/treemix/populations.snps.lFilt.iFilt.input.cov.gz"
#     output:
#         "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/treemix/populations.snps.lFilt.iFilt.input.tree.pdf"
#     params:
#         treemix_dir = config["treemixDir"],
#         m = config["treemixParams"]["m"],
#     resources:
#         mem_mb = 8000,
#         runtime = 720,
#         tmpdir = "results/tmp",
#         partition = "serial_requeue,shared"
#     shell:
#         """
#         module load gcc/12.2.0-fasrc01 
#         module load R/4.2.2-fasrc01

#         out_prefix={output}
#         out_prefix=${{out_prefix%.cov*}}
        
        
#         """
