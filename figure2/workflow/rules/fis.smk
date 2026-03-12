###################
# FIS #
###################

# Compute FIS values
rule run_fis:
    input:
        vcf = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.vcf.gz",
        popmap = "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/populations.all.lFilt.iFilt.pop"
    output:
        "results/{breadth}/{sp_grp}" + f"/res_stacks/{paramspace.wildcard_pattern}" + "/mind_{mind}/fis/populations.all.lFilt.iFilt_fis.txt"
    params:
        wsize = config["pixyParams"]["windowSize"],
        minsize = config["pixyParams"]["minChromSize"]
    resources:
        mem_mb_per_cpu = 10000,
        runtime = 480,
        tmpdir = "/scratch",
        slurm_partition = "shared"
    threads: 8
    conda:
        "../envs/pixy.yaml"
    shell:
        """
        touch {input.vcf}.tbi

        out_prefix={output}
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

        ./workflow/scripts/FIS.py \
            --vcf {input.vcf} \
            --populations $pop_fn \
            --window_size {params.wsize} \
            --n_cores {threads} \
            --chromosomes $chromstr \
            --outfile {output}
        """


