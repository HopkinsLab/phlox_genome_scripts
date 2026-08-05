import yaml
import glob, os, pathlib


#Config file
configfile: "config/config.yaml"
res_config = yaml.load(open("config/resources.yaml"),Loader=yaml.FullLoader)
samsheet="config/samples.tsv"

samdict = {}
with open(samsheet) as fin:
    for line in fin:
        line = line.rstrip()
        data = line.split("\t")
        samdict[data[0]] = {'bamfile' : data[1]}

samples = list(samdict.keys())

#Global variables from config
ref_org = config["ref_org"]
int_files = config["intermediate_files"]

#Functions
def get_bam(wildcards):
    bfile = samdict[wildcards.sample]['bamfile']
    return bfile


#Global rule
rule all:
    input:
        expand("results/{sample}_finalCall.sample.vcf", sample=samples),
        "results/JAS_MERGED_ALLSAMPLES.vcf"


#Rules
rule snifflesCall:
    input:
        ref = config["reference"],
        bamfile = get_bam
    output:
        vcf = int_files + "{sample}_vs_" + ref_org + "_sniffles.vcf",
        snf = int_files +  "{sample}_vs_" + ref_org + "_sniffles.snf"
    conda:
        "workflow/envs/sniffles.yaml"
    resources:
        mem_mb = res_config["sniffles_call"]["mem_mb"],
        time = res_config["sniffles_call"]["time"]
    threads:
        res_config['sniffles_call']['threads']
    shell:
        "sniffles --threads {threads} --reference {input.ref} --input {input.bamfile} --vcf {output.vcf} --snf {output.snf}"

rule cutesvCall:
    input:
        ref = config["reference"],
        bamfile = get_bam
    output:
        vcf = int_files + "{sample}_vs_" + ref_org + "_cutesv.vcf"
    conda:
        "workflow/envs/cutesv.yaml"
    resources:
        mem_mb = res_config["cutesv_call"]["mem_mb"],
        time = res_config["cutesv_call"]["time"]
    threads:
        res_config['cutesv_call']['threads']
    params:
        tmpdir = int_files + "{sample}"
    shell:
        """
        mkdir -p {params.tmpdir}
        cuteSV -t {threads} --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DE 0.3 {input.bamfile} {input.ref} {output.vcf} {params.tmpdir}
        """

rule svimCall:
    input:
        ref = config["reference"],
        bamfile = get_bam
    output:
        vcf = int_files + "{sample}_vs_" + ref_org + "_svim.vcf"
    conda:
        "workflow/envs/svim.yaml"
    resources:
        mem_mb = res_config["svim_call"]["mem_mb"],
        time = res_config["svim_call"]["time"]
    params:
        outdir = "workflow/{sample}_svim"
    shell:
        """
        svim alignment {params.outdir} {input.bamfile} {input.ref}
        bcftools view -i 'QUAL >= 10' {params.outdir}/variants.vcf > {output.vcf}
        """

rule survMergeCallers:
    input:
        snif = int_files + "{sample}_vs_" + ref_org + "_sniffles.vcf",
        cute = int_files + "{sample}_vs_" + ref_org + "_cutesv.vcf",
        svim = int_files + "{sample}_vs_" + ref_org + "_svim.vcf"
    output:
        int_files + "{sample}_merged.vcf"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1.5 * res_config["vcf_merge"]["mem_mb"],
        time = res_config["vcf_merge"]["time"]
    params:
        surv_path = config["survivor_path"],
        text = "{sample}_vcfs.txt"
    shell:
        """
        ls {input.snif} {input.cute} {input.svim} > {params.text}
        {params.surv_path} merge {params.text} 1000 3 1 1 0 50 {output}
        """

rule jasMergeCallers:
    input:
        snif = int_files + "{sample}_vs_" + ref_org + "_sniffles.vcf",
        cute = int_files + "{sample}_vs_" + ref_org + "_cutesv.vcf",
        svim = int_files + "{sample}_vs_" + ref_org + "_svim.vcf"
    output:
        int_files + "{sample}_jas_merge_callers.vcf"
    conda:
        "workflow/envs/jasmine.yaml"
    resources:
        mem_mb = res_config["jas_merge"]["mem_mb"],
        time = res_config["jas_merge"]["time"]
    params:
        text = "{sample}_vcfs2.txt"
    threads:
        res_config['jas_merge']['threads']
    shell:
        """
        ls {input.snif} {input.cute} {input.svim} > {params.text}
        jasmine file_list={params.text} out_file={output} --allow_intrasample threads={threads}
        """

rule vcfFilter:
    input:
        surv = int_files + "{sample}_merged.vcf",
        jas = int_files + "{sample}_jas_merge_callers.vcf"
    output:
        survfilt = int_files + "{sample}_merged_filtered.vcf",
        jasfilt = int_files + "{sample}_jas_merge_filtered.vcf"
    resources:
        mem_mb = res_config["vcf_filter"]["mem_mb"],
        time = res_config["vcf_filter"]["time"]
    params:
        gaps = config["gap_bed"],
        surv_path = config["survivor_path"],
 #       survinter = int_files + "{sample}_merged_filtered.vcf",
 #       jasinter = int_files + "{sample}_jas_merge_filtered.vcf"
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        """
        grep -v '<TRA>' {input.surv} | grep -v 'BND' | bcftools view -i 'SVLEN<100000' - > {output.survfilt}
        grep -v '<TRA>' {input.jas} | grep -v 'BND' | bcftools view -i 'SVLEN<100000' - > {output.jasfilt}
        """

 #       {params.surv_path} filter {params.survinter} {params.gaps} 50 -1 0 -1 {output.survfilt}
 #       {params.surv_path} filter {params.jasinter} {params.gaps} 50 -1 0 -1 {output.jasfilt}

rule survMergeSamples:
    input:
        vcfs = expand(int_files + "{sample}_merged_filtered.vcf", sample=samples),
    output:
        mergevcf = "results/ALL_MERGED.vcf"
    resources:
        mem_mb = res_config["vcf_merge"]["mem_mb"],
        time = res_config["vcf_merge"]["time"]
    params:
        surv_path = config["survivor_path"],
        text = "all_vcfs.txt"
    shell:
        """
        ls {input.vcfs} > {params.text}
        {params.surv_path} merge {params.text} 1000 1 1 1 0 50 {output.mergevcf}
        """

rule jasMergeSamples:
    input:
        vcfs = expand(int_files + "{sample}_jas_merge_filtered.vcf", sample=samples),
    output:
        "results/JAS_MERGED_ALLSAMPLES.vcf"
    conda:
        "workflow/envs/jasmine.yaml"
    resources:
        mem_mb = res_config["jas_merge"]["mem_mb"],
        time = res_config["jas_merge"]["time"]
    params:
        text = "all_vcfs2.txt"
    threads:
        res_config['jas_merge']['threads']
    shell:
        """
        ls {input.vcfs} > {params.text}
        jasmine file_list={params.text} out_file={output} --allow_intrasample --output_genotypes threads={threads}
        """

rule forceCall:
    input:
        ref = config["reference"],
        bamfile = get_bam,
        knownsv = "results/ALL_MERGED.vcf"
    output:
        vcf = "results/{sample}_finalCall.vcf"
    resources:
        mem_mb = res_config["force_call"]["mem_mb"],
        time = res_config["force_call"]["time"]
    threads:
        res_config['force_call']['threads']
    conda:
        "workflow/envs/sniffles.yaml"
    shell:
        "sniffles --threads {threads} --input {input.bamfile} --genotype-vcf {input.knownsv} --vcf {output.vcf}"


rule fixHeader:
    input:
        "results/{sample}_finalCall.vcf"
    output:
        vcf =  "results/{sample}_finalCall.sample.vcf"
    resources:
        mem_mb = res_config["vcf_filter"]["mem_mb"],
        time = res_config["vcf_filter"]["time"]
    conda:
        "workflow/envs/bcftools.yaml"
    params:
        rename = temp("{sample}_rename.txt"),
        trimvcf = temp("{sample}_trimmed.vcf")
    shell:
        """
        echo "{wildcards.sample}" > {params.rename}
        cut -f1-10 {input} > {params.trimvcf}
        bcftools reheader -s {params.rename} {params.trimvcf} -o {output.vcf}
        """

