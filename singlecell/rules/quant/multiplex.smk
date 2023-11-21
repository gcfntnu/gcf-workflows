
#from os.path import join

include:
    'doublets.smk'

#container: 'docker://drneavin/demuxafy:v2.0.1'
#container: 'Demuxafy.sif'
    
def get_barcodes(wildcards):
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', True):
            return rules.cellranger_cellbender.output.csv
        return rules.cellranger.output.filtered_barcodes
    elif wildcards.quantifier == 'starsolo':
        return rules.starsolo_quant.output.barcodes
    else:
        pass

def get_bam(wildcards):
    #join(QUANT_INTERIM, 'star', '{sample}', '{sample}_Aligned.sortedByCoord.out.bam')
    #join(QUANT_INTERIM, 'cellranger', '{sample}', 'outs', '{sample}_possorted_genome_bam.bam')
    if wildcards.quantifier == 'cellranger':
        return rules.cellranger_quant.output.bam
    elif wildcards.quantifier == 'starsolo':
        return rules.starsolo_bam.output
    else:
        pass

def get_common_genotypes_vcf():
    return 'data/ext/common_variants_grch38.sorted.vcf'

def get_donor_bam(wildcards):
    donor_files = config['samples'][wildcards.sample]['donor_files'].split(',')
    return [os.path.join(config['quant']['demultiplex'].get('donor_dir', ''), fn) for fn in donor_files]

def get_donor_ids(wildcards):
    donor_ids = config['samples'][wildcards.sample]['donor_ids'].split(',')
    return donor_ids

rule bam_rename_chromosomes:
    input:
        bam = get_bam,
        mapping_file = '/mnt/archive/ext_cache/ext/ensembl/release-109/homo_sapiens/GRCh38/anno/ensembl_to_10xgenomics_2020-A_map.txt'
    output:
        bam = '{quantifier}/{sample}_fixed.bam'
    container:
        'docker://lindenb/jvarkit:1b2aedf24'
    shell:
        'java -jar /opt/jvarkit/dist/jvarkit.jar bamrenamechr '
        '-f {input.mapping_file} '
        '--samoutputformat BAM '
        '--out {output.bam} '
        '{input.bam} '

rule freemuxlet_pileup:
    input:
        bam = get_bam,
        barcodes = get_barcodes,
        vcf = get_common_genotypes_vcf(),
    output:
        cel = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'freemuxlet', 'pileup.cel.gz')
    params:
        freemuxlet_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'freemuxlet')
    threads:
        24
    container:
        'docker://yenchungchen/popscle:latest'
    shell:
        'popscle dsc-pileup '
        '--sam {input.bam} '
        '--vcf {input.vcf} '
        '--group-list {input.barcodes} '
        '--out {params.freemuxlet_dir}/pileup '
        
rule cellsnp_pileup_1b: # pileup with defined snps (bulk)
    input:
        bam = get_donor_bam,
        vcf = get_common_genotypes_vcf()
    output:
        vcf_base = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.base.vcf.gz'),
        vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.vcf.gz'),
        samples = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.samples.tsv'),
        mtx_ad = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.tag.AD.mtx'),
        mtx_dp = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.tag.DP.mtx'),
        mtx_other = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.tag.OTH.mtx')
    params:
        cellsnp_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp'),
        input_bam = lambda wildcards, input: ','.join(input.bam),
        donor_ids = lambda wildcards: ','.join(get_donor_ids(wildcards)),
        p = 20,
        minMAF = 0.1,
        minCOUNT = 20,
        minLEN = 30,
        minMAPQ = 20,
        args = '--gzip --genotype --cellTAG None --UMItag None '
    threads:
        16
    container:
        'docker://gcfntnu/cellsnp-lite:1.2.3'
    shell:
        'cellsnp-lite '
        '-s {params.input_bam} '
        '-I {params.donor_ids} '
        '-R {input.vcf} '
        '-O {params.cellsnp_dir} '
        '-p {params.p} '
        '--minMAF {params.minMAF} '
        '--minCOUNT {params.minCOUNT} ' 
        ' {params.args} '
        '--nproc {threads} '

rule index_donor_vcf:
    input:
        donor_vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.vcf.gz')
    output:
        donor_vcf_index = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.vcf.gz.tbi')
    singularity:
        'docker://gcfntnu/vcftools:1.18'
    shell:
        'tabix {input} '
        
if config['db']['reference_db'] == '10xgenomics':
    rule cellsnp_donor_vcf:
        input:
            vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.vcf.gz'),
            vcf_index = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.vcf.gz.tbi'),
            chrom_map = '/mnt/archive/ext_cache/ext/ensembl/release-109/homo_sapiens/GRCh38/anno/ensembl_to_10xgenomics_2020-A_map.txt'
        output:
            vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.chr.vcf.gz')
        singularity:
            'docker://gcfntnu/vcftools:1.18'
        shell:
            'bcftools annotate --rename-chrs {input.chrom_map} {input.vcf} -Oz -o {output.vcf} '
else:
    rule cellsnp_donor_vcf:
        input:
            vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.vcf.gz'),
            vcf_index = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.vcf.gz.tbi'),
            chrom_map = '/mnt/archive/ext_cache/ext/ensembl/release-109/homo_sapiens/GRCh38/anno/ensembl_to_10xgenomics_2020-A_map.txt'
        output:
            vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.chr.vcf.gz')
        singularity:
            'docker://gcfntnu/vcftools:1.18'
        shell:
            'ln -s {input.vcf} {output.vcf} && ln -s {input.vcf_index} {output.vcf_index}'

rule cellsnp_pileup_1a: # pileup with defined snps and barcodes (singlecell)
    input:
        bam = get_bam,
        barcodes = get_barcodes,
        vcf = get_common_genotypes_vcf()
    output:
        vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'cellsnp', 'cellSNP.base.vcf.gz'),
        samples = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'cellsnp', 'cellSNP.samples.tsv'),
        mtx_ad = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'cellsnp', 'cellSNP.tag.AD.mtx'),
        mtx_dp = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'cellsnp', 'cellSNP.tag.DP.mtx'),
        mtx_other = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'cellsnp', 'cellSNP.tag.OTH.mtx')
    params:
        cellsnp_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'cellsnp'),
        input_bam = lambda wildcards, input: input.bam if isinstance(input.bam, str) else ','.join(input.bam),
        p = 20,
        minMAF = 0.1,
        minCOUNT = 20,
        minLEN = 30,
        minMAPQ = 20,
        args = '--gzip '
    threads:
        16
    shell:
        'cellsnp-lite '
        
        '-s {params.input_bam} '
        '-b {input.barcodes} '
        '-R {input.vcf} '
        '-O {params.cellsnp_dir} '
        '-p {params.p} '
        '--minMAF {params.minMAF} '
        '--minCOUNT {params.minCOUNT} ' 
        ' {params.args} '
        '--nproc {threads} '

rule vireo_noref:
    input:
        rules.cellsnp_pileup_1a.output
    output:
        donor_ids = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'vireo', 'donor_ids.tsv'),
        doublet = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'vireo', 'prob_doublet.tsv.gz'),
        singlet = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'vireo', 'prob_singlet.tsv.gz'),
        summary = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'vireo', 'summary.tsv')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        vireo_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'vireo'),
        cellsnp_dir= rules.cellsnp_pileup_1a.params.cellsnp_dir
    shell:
        'vireo '
        '-c {params.cellsnp_dir} '
        '-o {params.vireo_dir} '
        '-N {params.n} '
    
rule vireo_ref:
    input:
        donor_vcf = rules.cellsnp_donor_vcf.output.vcf,
        single_cell_vcf = rules.cellsnp_pileup_1a.output
    output:
        donor_ids = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing', 'with_reference', 'vireo', 'donor_ids.tsv'),
        doublet = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing', 'with_reference', 'vireo', 'prob_doublet.tsv.gz'),
        singlet = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing', 'with_reference', 'vireo', 'prob_singlet.tsv.gz'),
        summary = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing', 'with_reference', 'vireo', 'summary.tsv')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        vireo_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'with_reference', 'vireo'),
        cellsnp_dir= rules.cellsnp_pileup_1a.params.cellsnp_dir
    threads:
        12
    shell:
        'vireo '
        '-d {input.donor_vcf} '
        '-c {params.cellsnp_dir} '
        '-o {params.vireo_dir} '
        '-N {params.n} '
        '--callAmbientRNAs '
        '--nproc {threads} '

rule souporcell_noref:
    input:
        bam = get_bam,
        barcodes = get_barcodes,
        vcf = get_common_genotypes_vcf(),
        genome = join(REF_DIR, 'fasta', 'genome.fa')
    output:
        vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'souporcell', 'cluster_genotypes.vcf')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        souporcell_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'souporcell'),
        args = '--ploidy 2 --min_alt 10 --min_ref 10 --max_loci 2048 --skip_remap True '
    threads:
        24
    container:
        'shub://wheaton5/souporcell:latest'
    shell:
        'souporcell_pipeline.py '
        '-i {input.bam} '
        '-b {input.barcodes} '
        '-f {input.genome} '
        '--common_variants {input.vcf} '
        '-t {threads} '
        '-o {params.souporcell_dir} '
        '-k {params.n} '
        '{params.args} '


rule decompressed_vcf:
    input:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.chr.vcf.gz')
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.chr.vcf')
    container:
        'docker://gcfntnu/vcftools:1.18'
    shell:
        'bgzip -d -c {input} > {output}'


rule souporcell_ref:
    input:
        bam = get_bam,
        barcodes = get_barcodes,
        donor_vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'donor_cellsnp', 'cellSNP.cells.chr.vcf'),
        genome = join(REF_DIR, 'fasta', 'genome.fa'),
    output:
        vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'with_reference', 'souporcell', 'cluster_genotypes.vcf')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        souporcell_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'with_reference', 'souporcell'),
        donor_ids = lambda wildcards: ' '.join(get_donor_ids(wildcards)),
        args = '--skip_remap True '
    threads:
        24
    container:
        'shub://wheaton5/souporcell:latest'
    shell:
        'souporcell_pipeline.py '
        '-i {input.bam} '
        '-b {input.barcodes} '
        '-f {input.genome} '
        '--known_genotypes {input.donor_vcf} '
        '-t {threads} '
        '-o {params.souporcell_dir} '
        '-k {params.n} '
        '{params.args} '
        '--known_genotypes_sample_names {params.donor_ids} '

rule freemuxlet_noref:
    input:
        barcodes = get_barcodes,
        cel = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'freemuxlet', 'pileup.cel.gz')
    output:
        vcf = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'freemuxlet', 'freemuxlet.clust1.vcf.gz')
    params:
        freemuxlet_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'freemuxlet'),
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4)
    threads:
        24
    container:
        'docker://yenchungchen/popscle:latest'
    shell:
        'popscle freemuxlet '
        '--group-list {input.barcodes} '
        '--plp {params.freemuxlet_dir}/pileup '
        '--out {params.freemuxlet_dir}/freemuxlet '
        '--nsample {params.n}'

rule demuxlet_ref:
    input:
        barcodes = get_barcodes,
        cel = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'freemuxlet', 'pileup.cel.gz'),
        donor_vcf = rules.cellsnp_donor_vcf.output.vcf
    output:
        join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'with_reference', 'demuxlet', 'demuxlet.best')
    params:
        freemuxlet_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'freemuxlet'),
        demuxlet_prefix = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'with_reference', 'demuxlet', 'demuxlet'),
        demuxlet_dir = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing',  'with_reference', 'demuxlet'),
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        field = 'GT',
        args = ' '
    threads:
        24
    container:
        'docker://yenchungchen/popscle:latest'
    shell:
        'popscle demuxlet '
        '--plp {params.freemuxlet_dir}/pileup '
        '--vcf {input.donor_vcf} '
        '--field {params.field} '
        '--group-list {input.barcodes} '
        '--out {params.demuxlet_prefix} '
        '{params.args}'



def get_demuxafy_methods(test=False):
    """Preprint suggestions.

    needs update when I can get souporcell with reference to work
    """
    if config['quant'].get('demultiplex'):
        n = config['quant']['demultiplex'].get('n', 4)
        if config['quant']['demultiplex']['with_reference']:
            if n < 4:
                methods = ['demuxlet', 'vireo_ref', 'scds']
            elif n < 16:
                methods = ['demuxlet', 'vireo_ref', 'scds']
            elif n < 64:
                methods = ['demuxlet', 'vireo_ref', 'scds']
            else:
                methods = ['demuxlet', 'vireo_ref', 'scds']
        else:
            if n < 4:
                methods = ['freemuxlet', 'souporcell', 'scds']
            elif n < 16:
                methods = ['freemuxlet', 'souporcell', 'vireo']
            elif n < 64:
                methods = ['freemuxlet', 'souporcell', 'vireo']
            else:
                methods = ['freemuxlet', 'souporcell']
    else:
        methods = ['solo', 'doublet_detection', 'scds']
    return methods


def get_demuxafy_method_args():
    args = ''
    for m in get_demuxafy_methods():
        if m == 'freemuxlet':
            args += ' --freemuxlet ' + rules.freemuxlet_noref.params.freemuxlet_dir
        elif m == 'souporcell':
            args += ' --souporcell ' + rules.souporcell_noref.params.souporcell_dir
        elif m == 'vireo':
            args += ' --vireo ' + rules.vireo_noref.params.vireo_dir
        elif m == 'vireo_ref':
            args += ' --vireo ' + rules.vireo_ref.params.vireo_dir
        elif m == 'demuxlet':
            args += ' --demuxlet ' + rules.demuxlet_ref.params.demuxlet_dir
        elif m == 'scds':
            args += ' --scds ' + rules.Scds.params.out_dir
        elif m == 'solo':
            args += ' --solo ' + rules.solo_summary.params.out_dir
        elif m == 'doublet_detection':
            args += ' --DoubletDetection ' + rules.DoubletDetection.params.out_dir
        elif m == 'scdbl_finder':
            args += ' --scDblFinder ' + rules.ScDblFinder.params.out_dir
        elif m == 'scrublet':
            args += ' --scrublet ' + rules.scrublet.params.out_dir
        
    return args

def get_aggr_input(wildcards):
    input = {}
    for m in get_demuxafy_methods():
        if m == 'fremuxlet':
            input['freemuxlet'] = rules.freemuxlet_noref.output[0]
        elif m == 'demuxlet':
            input['demuxlet'] = rules.demuxlet_ref.output[0]
        elif m == 'souporcell':
            input['souporcell'] = rules.souporcell_noref.output[0]
        elif m == 'vireo':
            input['vireo'] = rules.vireo_noref.output[0]
        elif m == 'vireo_ref':
            input['vireo_ref'] = rules.vireo_ref.output[0]
        elif m == 'scds':
            input['scds'] = rules.Scds.output[0]
        elif m == 'solo':
            input['solo'] = rules.solo_summary.output[0]
        elif m == 'doublet_detection':
            input['doublet_detection'] = rules.DoubletDetection.output[0]
        elif m == 'scdbl_finder':
            input['scdbl_finder'] = rules.ScDblFinder.output[0]
        elif m == 'scrublet':
            input['scrublet'] = rules.scrublet.output[0]
    return input      

rule sample_aggregate_demux_results:
    input:
        unpack(get_aggr_input)
    params:
        method_dirs = get_demuxafy_method_args(),
        combination = "MajoritySinglet",
        ref_arg = '--ref vireo' if 'vireo' in get_demuxafy_methods() else ''
    container:
        'Demuxafy.sif'
    output:
        combined = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing', 'with_reference', 'combined_demux.txt'),
        combined_assigned = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing', 'with_reference', 'combined_demux_w_combined_assignments.tsv'),
        upset_fig = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing', 'with_reference', 'combined_demuxSinglets_upset.pdf')
    shell:
        'Combine_Results.R '
        '{params.method_dirs} '
        '--method {params.combination} '
        '--out {output.combined} '



demultiplexer = config['quant'].get('demultiplex', {}).get('method')
quantifier = config['quant']['method']
if demultiplexer == 'demuxafy':
    rule multiplex_all:
        input:
            expand(rules.sample_aggregate_demux_results.output.combined_assigned, quantifier=config['quant']['method'], sample=SAMPLES)
        output:
            join(QUANT_INTERIM, 'aggregate',config['quant']['method'] , 'droplet_type.txt')
        params:
            script = srcdir("scripts/combine_demultiplex.py")
        container:
            'docker://' + config['docker']['default']
        shell:
            'python {params.script} '
            '{input} '
            '-o {output} '

elif demultiplexer == 'vireo_ref':
    rule multiplex_all:
        input:
             expand(rules.vireo_ref.output, sample=SAMPLES, quantifier=quantifier)
