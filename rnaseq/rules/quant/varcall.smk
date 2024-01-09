# -*- mode: snakemake -*-
"""

fixme: cellsnp-lite lacks implmentation of non-aggregated path
fixme: cellsnp-lite assumes reference_db is ensembl (or has ensemble type chromosome names with no chr-prefix)
"""
VARCALL_INTERIM = join(ALIGN_INTERIM, ALIGNER, 'varcall')

rule gene_intervals:
    input:
        bed = join(REF_DIR, 'anno', 'genes.bed'),
        genome = join(REF_DIR, 'fasta', 'genome.dict')
    output:
        join(REF_DIR, 'anno', 'genes.intervals')
    singularity:
        'docker://' + config['docker']['picard_gatk']
    shell:
        'gatk BedToIntervalList '
        '-I={input.genome} '
        '-O={output} '
        '-SD={input.genome} '

rule picard_add_rg:
    input:
        bam = get_sorted_bam
    params:
        machine = config.get('machine', 'nextseq500').replace(' ', '')
    output:
        bam = join(VARCALL_INTERIM, '{sample}_RG.bam')
    singularity:
        'docker://' + config['docker']['picard_gatk']
    shell:
        'picard AddOrReplaceReadGroups '
        'INPUT={input.bam} '
        'OUTPUT={output.bam} '
        'RGID={wildcards.sample} '
        'RGLB={wildcards.sample} '
        'RGPL=illumina '
        'RGPU={params.machine} '
        'RGSM={wildcards.sample} '
    
rule gatk_split_and_trim:
    input:
        bam = join(VARCALL_INTERIM, '{sample}_RG.bam'),
        genome = join(REF_DIR, 'fasta', 'genome.fa')
    output:
        bam = join(VARCALL_INTERIM, '{sample}.bam.split')
    params:
        '-RF ReassignOneMappingQuality -RMQF 255 -RMQT 60 '
    singularity:
        'docker://' + config['docker']['picard_gatk']
    threads:
        4
    shell:
        'gatk SplitNCigarReads '
        '-R {input.genome} '
        '-I {input.bam} '
        '-O {output.bam} '

rule gatk_varcall:
    input:
        bam = join(VARCALL_INTERIM, '{sample}.bam.split'),
        genome = join(REF_DIR, 'fasta', 'genome.fa')
    output:
        vcf = join(VARCALL_INTERIM, '{sample}_all.vcf'),
        bam = join(VARCALL_INTERIM, '{sample}.bam')
    params:
        '-stand-call-conf 10.0 '
    singularity:
        'docker://' + config['docker']['picard_gatk']
    threads:
        48
    shell:
        'gatk HaplotypeCaller '
        '--native-pair-hmm-threads {threads} '
        '-R {input.genome} '
        '-I {input.bam} '
        '-O {output.vcf} '
        '-bamout {output.bam} '

rule gatk_var_filter:
    input:
        vcf = join(VARCALL_INTERIM, '{sample}_all.vcf'),
        genome = join(REF_DIR, 'fasta', 'genome.fa')
    output:
        vcf = join(VARCALL_INTERIM, '{sample}.filtered.vcf')
    singularity:
        'docker://' + config['docker']['picard_gatk']
    shell:
        """
        gatk VariantFiltration --R {input.genome} --V {input.vcf} --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" -O {output.vcf}
        """

rule gatk_varcall_all:
    input:
        expand(join(VARCALL_INTERIM, '{sample}.filtered.vcf'), sample=SAMPLES)
    output:
        vcf = join(VARCALL_INTERIM, 'merged.vcf')
    params:
        vcfs = lambda wildcards, input : '--INPUT ' + ' --INPUT  '.join(input)
    singularity:
        'docker://' + config['docker']['picard_gatk']
    shell:
        'gatk MergeVcfs '
        '{params.vcfs} '
        '--OUTPUT {output} '

rule donor_list:
    output:
        join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'donor_list.txt')
    params:
        donor_ids = lambda wildcards: r'\n'.join(AGGR_IDS.get(wildcards.aggr_id, []))
    shell:
        'echo -e "{params.donor_ids}" > {output}'

def get_aggr_bam(wildcards):
    sub_samples = AGGR_IDS[wildcards.aggr_id]
    ALIGNER = config.get('align', {}).get('genome', {}).get('aligner', 'star')
    return expand(join(ALIGN_INTERIM, ALIGNER, '{sample}.sorted.bam'), sample=sub_samples)
    
rule cellsnp_pileup_1b_aggr: # pileup with defined variants (bulk)
    input:
        bam = get_aggr_bam,
        vcf = join(REF_DIR, 'anno', 'common_variants.vcf'),
        donor_list = join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'donor_list.txt')
    output:
        vcf_base = join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.base.vcf'),
        vcf = join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.vcf'),
        samples = join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.samples.tsv'),
        mtx_ad = join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.tag.AD.mtx'),
        mtx_dp = join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.tag.DP.mtx'),
        mtx_other = join(VARCALL_INTERIM, 'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.tag.OTH.mtx')
    params:
        cellsnp_dir = join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}'),
        input_bam = lambda wildcards, input: ','.join(input.bam),
        minMAF = 0.1,
        minCOUNT = 20,
        args = '--genotype --cellTAG None --UMItag None '
    threads:
        24
    container:
        'docker://' + config['docker']['cellsnp-lite']
    shell:
        'cellsnp-lite '
        '-s {params.input_bam} '
        '-i {input.donor_list} '
        '-R {input.vcf} '
        '-O {params.cellsnp_dir} '
        '--minMAF {params.minMAF} '
        '--minCOUNT {params.minCOUNT} '
        ' {params.args} '
        '--nproc {threads} '

rule index_donor_vcf:
    input:
        donor_vcf = join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.vcf')
    output:
        donor_vcf_index = join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.vcf.gz.tbi'),
        donor_vcf =  join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.vcf.gz')
    container:
        'docker://' + config['docker']['vcftools']
    threads:
        4
    shell:
        'bgzip --threads {threads} --stdout --index --index-name {output.donor_vcf_index} {input.donor_vcf} > {output.donor_vcf} '

rule chrom_conversion_chr:
    output:
        txt = join(REF_DIR, 'anno', 'chrom_map.txt')
    run:
        chrom = [str(i) for i in range(1, 23)]
        chrom.extend(['X', 'Y', 'MT'])
        with open(output.txt, 'w') as fh:
            for c in chrom:
                line = '{}\t{}\n'.format(c, 'chr'+c)
                fh.write(line)

rule cellsnp_donor_vcf_chr:
    input:
        vcf = join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.vcf.gz'),
        chrom_map = join(REF_DIR, 'anno', 'chrom_map.txt')
    output:
        vcf = join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.chr.vcf')
    singularity:
        'docker://' + config['docker']['vcftools']
    shell:
        'bcftools annotate --rename-chrs {input.chrom_map} {input.vcf} -Ov -o {output.vcf} '
                                      
rule index_donor_vcf_chr:
    input:
        donor_vcf = join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.chr.vcf')
    output:
        donor_vcf_index = join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.chr.vcf.gz.tbi'),
        donor_vcf =  join(VARCALL_INTERIM,  'aggregate', 'cellsnp', '{aggr_id}', 'cellSNP.cells.chr.vcf.gz')
    container:
        'docker://' + config['docker']['vcftools']
    threads:
        4
    shell:
        'bgzip --threads {threads} --stdout --index --index-name {output.donor_vcf_index} {input.donor_vcf} > {output.donor_vcf} '

rule cellsnp_all:
        input:
            expand(rules.index_donor_vcf_chr.output.donor_vcf, aggr_id = AGGR_IDS)
                                      
                                      
