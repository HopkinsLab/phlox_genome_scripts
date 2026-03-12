###################
# PIXY #
###################

# Compute pairwise fst using pixy
rule run_pixy:
    input:
        vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.vcf.gz",
        popmap = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.pop"
    output:
        dxy = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pixy/populations.all.lFilt.iFilt_dxy.txt",
        pi = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pixy/populations.all.lFilt.iFilt_pi.txt",
        fst = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/pixy/populations.all.lFilt.iFilt_fst.txt"
    params:
        wsize = config["pixyParams"]["windowSize"],
        minsize = config["pixyParams"]["minChromSize"]
    resources:
        mem_mb_per_cpu = 3000,
        runtime = 240,
        tmpdir = "/scratch",
        slurm_partition = "serial_requeue"
    threads: 8
    conda:
        "../envs/pixy.yaml"
    shell:
        """
        touch {input.vcf}.tbi

        out_prefix={output.dxy}
        out_prefix=${{out_prefix%_dxy.txt}}

        pop_fn=$out_prefix.pop
        cp -f {input.popmap} $pop_fn

        # Chromosome string
        tmpfn=$(mktemp)
        tabix -H {input.vcf} | \
            grep '^##contig=' | \
            sed -e 's/^##contig=<ID=//' -e 's/,length=/ /' -e 's/>$//' | \
            awk -v ORS="," -v minsize={params.minsize} '$2 > minsize {{print $1}}' \
            > $tmpfn

        chromstr=$(cat $tmpfn)
        chromstr=${{chromstr%,}}

        pixy --stats pi fst dxy --vcf {input.vcf} --populations $pop_fn \
            --window_size {params.wsize} \
            --n_cores {threads} \
            --chromosomes $chromstr \
            --output_folder $( dirname  $out_prefix ) \
            --output_prefix $( basename $out_prefix )
        """


