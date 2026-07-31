localrules: make_gap_bed, prep_shuf_annot_intervals

# Download reference, make fasta index if needed, 
# simplify file naming
rule prep_reference:
    input:
        config["ref"]["fasta"]
    output:
        ref = "results/ref/" + config["ref"]["name"] + ".fasta",
        index = "results/ref/" + config["ref"]["name"] + ".fasta.fai",
    resources:
        mem_mb = 4000,
        runtime = 120,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        cp {input} {output.ref}
        if [[ -s {input}.fai ]]; then
            cp {input}.fai {output.index}
        else
            samtools faidx {output.ref}
        fi
        """

# Make windows 
rule make_windows:
    input:
        "results/ref/" + config["ref"]["name"] + ".fasta.fai"
    output:
        "results/ref/" + config["ref"]["name"] + ".w_{wname}.bed"
    # params:
    #     window_size = config["windowSize"]
    resources:
        mem_mb = 4000,
        runtime = 120,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        module load R/4.2.2-fasrc01; module load gcc/12.2.0-fasrc01
        wsize=$( workflow/scripts/wsize2wname.R {wildcards.wname} )
        bedtools makewindows -g {input} -w $wsize > {output}
        """

# convert gff annotations to bed
def get_gff_path(wildcards):
    return config["gffs"][wildcards.annot_type]["path"]


rule gffs2bed:
    input:
        fai = "results/ref/" + config["ref"]["name"] + ".fasta.fai",
        gff = get_gff_path
    output:
        "results/annotations/{annot_type}.bed"
    params:
        filter_command = lambda wildcards: config["gffs"][wildcards.annot_type]["filterCommand"]
    resources:
        mem_mb = 8000,
        runtime = 120,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        ext={input.gff}
        ext=${{ext##*.}}
        if [[ $ext == "gz" ]]; then
            catcmd="zcat"
        else
            catcmd="cat"
        fi

        $catcmd {input.gff} | \
            {params.filter_command} | \
            gff2bed --do-not-sort | \
            bedtools sort -i /dev/stdin -g {input.fai} \
            > {output}
        """

rule annot_size:
    input:
        gaps = "results/annotations/gaps.bed",
        annot = "results/annotations/{annot_type}.bed"
    output:
        "results/annotations/{annot_type}.uniq_bp.txt"
    resources:
        mem_mb = 2000,
        runtime = 20,
        tmpdir = "/scratch",
        partition = "serial_requeue"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        tmpfn=$(mktemp --suffix=.bed)
        sort -k1,1 -k2,2n {input.gaps} > $tmpfn

        sort -k1,1 -k2,2n {input.annot} | \
            bedtools merge -i - | \
            bedtools subtract -a - -b $tmpfn | \
            awk '{{sum += $3 - $2}} END {{print sum}}' \
            > {output}
        """


# annotate repeat
rule annotate_repeat_families:
    input:
        bed = "results/annotations/repeats.bed",
        rpt_fai = config["ref"]["rpt_fai"]
    output:
        "results/annotations/repeats.families.bed"
    resources:
        mem_mb = 4000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    shell:
        """
        module load R/4.2.2-fasrc01; module load gcc/12.2.0-fasrc01
        workflow/scripts/annotate_repeat_families.R {input.bed} {input.rpt_fai}
        """

rule merge_repeats:
    input:
        "results/annotations/repeats.families.bed"
    output:
        "results/annotations/repeats.families.merged.bed"
    resources:
        mem_mb = 2000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        sort -k1,1 -k2,2n {input} | \
            bedtools merge -i - > {output}
        """


# Create bed file of gap locations
rule make_gap_bed:
    input:
        config["ref"]["gap_bed"]
    output:
        "results/annotations/gaps.bed"
    resources:
        mem_mb = 2000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    shell:
        """
        cp {input} {output}
        """

rule make_svs_bed:
    input:
        config["svs"]["path"]
    output:
        "results/annotations/SV_{sv_type}.{sv_sp}.bed"
    params:
        min_len = config["svs"]["minLen"],
        samples  = lambda wildcards: config["svs"]["samples"][wildcards.sv_sp]
    conda:
        "../envs/bedtools.yml"
    resources:
        mem_mb = 4000,
        runtime = 120,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    shell:
        """
        samp_str=""
        for s in {params.samples}; do
            samp_str="$samp_str$s,"
        done
        samp_str=${{samp_str%,}}

        bcftools view -s $samp_str {input} | \
            bcftools +fill-tags -- -t AN,AC | \
            bcftools query -i 'AC > 0' \
                -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\t%AC\n' | \
            awk -v OFS="\t" -v MINLEN={params.min_len} -v SVTYPE={wildcards.sv_type} \
                '$4 == SVTYPE {{if((($3 >= 0) && ($3 > MINLEN)) || (($3 < 0) && (-$3 > MINLEN))){{print $1,$2-1,$2,$5}}}}' \
            > {output}
        """

rule make_svs_qry:
    input:
        config["svs"]["path"]
    output:
        "results/annotations/SV.query.tsv"
    conda:
        "../envs/bedtools.yml"
    resources:
        mem_mb = 2000,
        runtime = 20,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    shell:
        """
        bcftools query \
            -Hf '%CHROM\t%POS\t%ID\t%SVLEN\t%SVTYPE[\t%GT]\n' {input} \
            > {output}
        """

rule sv_overlap:
    input:
        feat_bed = "results/annotations/{annot_type}.bed",
        sv_qry = "results/annotations/SV.query.tsv"
    output:
        "results/sv_overlap/{annot_type}.overlaps.bed",
    resources:
        mem_mb = 2000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        tmpbed=$(mktemp --suffix=.bed)
        grep -v '^#' {input.sv_qry} | \
            awk -v OFS="\t" '{{print $1,$2-1,$2,$3,$4,$5}}' > $tmpbed

        sort -k 1,1 -k2,2n {input.feat_bed} | \
            bedtools intersect -a $tmpbed -b - -u > {output}
        """


# coverage and counts across windows
rule windowed_coverage:
    input:
        ref_bed = "results/ref/" + config["ref"]["name"] + ".w_{wname}.bed",
        map_bed = "results/annotations/{annot_type}.bed"
    output:
        "results/coverage/{annot_type}.w_{wname}.bed"
    resources:
        mem_mb = 4000,
        runtime = 120,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedmap --delim '\t' --echo --bases-uniq --count \
            {input.ref_bed} {input.map_bed} | \
            awk -v OFS="\t" '{{print $1,$2,$3,".",".",".",$4,$5}}' > {output}
        """

# coverage and counts across windows
rule windowed_coverage_disjoint:
    input:
        ref_bed = "results/ref/" + config["ref"]["name"] + ".w_{wname}.bed",
        rpt_bed = "results/annotations/repeats.families.bed",
        exon_bed = "results/annotations/exons.bed",
        intron_bed = "results/annotations/introns.bed",
        gap_bed = "results/annotations/gaps.bed"
    output:
        touch("results/coverage/disjoint.w_{wname}.done")
    resources:
        mem_mb = 4000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "test" 
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        tmppre=$(mktemp)

        sort -k 1,1 -k2,2n {input.gap_bed} > $tmppre.gaps.bed
        sort -k 1,1 -k2,2n {input.rpt_bed} > $tmppre.rpt.bed

        sort -k 1,1 -k2,2n {input.exon_bed} | \
            bedtools subtract -a - -b $tmppre.gaps.bed | \
            sort -k 1,1 -k2,2n | \
            bedtools merge -i - \
            > $tmppre.exon.bed
        
        sort -k 1,1 -k2,2n {input.intron_bed} | \
            bedtools subtract -a - -b $tmppre.gaps.bed | \
            sort -k 1,1 -k2,2n | \
            bedtools merge -i - \
            > $tmppre.intron.bed

        awk -v FS="\t" -v OFS="\t" '$12 == "LTR"' $tmppre.rpt.bed | \
            bedtools subtract -a - -b $tmppre.exon.bed | \
            bedtools subtract -a - -b $tmppre.intron.bed | \
            bedtools subtract -a - -b $tmppre.gaps.bed | \
            sort -k 1,1 -k2,2n | \
            bedtools merge -i - \
            > $tmppre.repeats.LTR.bed

        awk -v FS="\t" -v OFS="\t" '$12 != "LTR"' $tmppre.rpt.bed | \
            bedtools subtract -a - -b $tmppre.exon.bed | \
            bedtools subtract -a - -b $tmppre.intron.bed | \
            bedtools subtract -a - -b $tmppre.gaps.bed | \
            bedtools subtract -a - -b $tmppre.repeats.LTR.bed | \
            sort -k 1,1 -k2,2n | \
            bedtools merge -i - \
            > $tmppre.repeats.other.bed


        tmp_ref=$(mktemp --suffix=".bed")
        sort-bed {input.ref_bed} > $tmp_ref

        tmp_fn=$(mktemp --suffix=".bed")


        outpre={output}
        outpre=${{outpre%.*}}

        nchar=${{#tmppre}}
        for fn in $(ls -1 $tmppre.*.bed | fgrep -v ".rpt.bed"); do
            suff=${{fn:$nchar}}
            sort-bed $fn > $tmp_fn
            bedmap --delim '\t' --echo --bases-uniq --count \
                        $tmp_ref $tmp_fn | \
                        awk -v OFS="\t" '{{print $1,$2,$3,".",".",".",$4,$5}}' \
                        > ${{outpre}}${{suff}}
        done
        """

######

checkpoint prep_shuf_annot_intervals:
    input:
        rpt_bed = "results/annotations/{annot_type}.bed",
    output:
        out_dir = directory("results/annot_shuffle/placeholders/{annot_type}/"),
        donefile = touch("results/annot_shuffle/placeholders/{annot_type}.done")
    params:
        shuf_reps = config["shufReps"]
    shell:
        """
        if [ ! -d {output.out_dir} ]; then
            mkdir -p {output.out_dir}
        fi

        for i in $(seq {params.shuf_reps}); do
            touch {output.out_dir}/$i.placeholder
        done
        """

rule shuf_annot_intervals:
    input:
        donefile = "results/annot_shuffle/placeholders/{annot_type}.done",
        annot_bed = "results/annotations/{annot_type}.bed",
        fai = "results/ref/" + config["ref"]["name"] + ".fasta.fai",
        scatter_placeholder = "results/annot_shuffle/placeholders/{annot_type}/{shuf_i}.placeholder"
    output:
        "results/annot_shuffle/bedfiles/{annot_type}/{shuf_i}.bed"
    resources:
        mem_mb = 2000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools shuffle -i {input.annot_bed} -g {input.fai} | sort -k1,1 -k2,2n > {output}
        """

rule shuf_rpt_overlap:
    input:
        shuf_bed = "results/annot_shuffle/bedfiles/{annot_type}/{shuf_i}.bed",
        rpt_bed = "results/annotations/repeats.families.merged.bed"
    output:
        bed = "results/annot_shuffle/bedfiles/{annot_type}/{shuf_i}.rpt_overlap.bed",
        summary = "results/annot_shuffle/bedfiles/{annot_type}/{shuf_i}.rpt_overlap.nbases.txt"
    resources:
        mem_mb = 2000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.shuf_bed} -b {input.rpt_bed} | \
            sort -k1,1 -k2,2n > {output.bed}

        bedtools merge -i {output.bed} | \
            awk '{{sum += $3 - $2}} END {{print sum}}' > {output.summary}
        """

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.prep_shuf_annot_intervals.get(annot_type=wildcards.annot_type).output[0]
    return expand("results/annot_shuffle/bedfiles/{annot_type}/{shuf_i}.rpt_overlap.nbases.txt",
                annot_type=wildcards.annot_type,
                shuf_i=glob_wildcards(os.path.join(checkpoint_output, "{shuf_i}.placeholder")).shuf_i)


rule pval_rpt_overlap:
    input:
        shuf_txts = aggregate_input,
        rpt_bed = "results/annotations/repeats.families.merged.bed",
        feat_bed = "results/annotations/{annot_type}.bed"
    output:
        bed = "results/annotations/{annot_type}.rpt_overlap.bed",
        summary = "results/annotations/{annot_type}.rpt_overlap.nbases.txt",
        pval = "results/annotations/{annot_type}.rpt_overlap.pval.txt",
    resources:
        mem_mb = 2000,
        runtime = 30,
        tmpdir = "/scratch",
        partition = "serial_requeue" 
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        sort -k 1,1 -k2,2n {input.feat_bed} | \
            bedtools intersect -a {input.rpt_bed} -b - | \
            sort -k 1,1 -k2,2n > {output.bed}

        bedtools merge -i {output.bed} | \
            awk '{{sum += $3 - $2}} END {{print sum}}' > {output.summary}

        nbases=$(cat {output.summary})
        cat {input.shuf_txts} | awk -v nbases=$nbases '{{sum += ($1 <= nbases)}} END {{print sum / NR}}' > {output.pval}
        """
