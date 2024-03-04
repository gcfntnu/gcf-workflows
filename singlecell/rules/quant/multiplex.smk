"""demultiplexing by natural variation

note: donor vcf files are assumed to contain the specifc samples to be demultiplexed and no more. 
"""

DEMUX_DIR = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing')


def get_singlecell_barcodes(wildcards):
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', False):
            return rules.cellranger_cellbender.output.aggr
        return rules.cellranger_quant.output.filt_barcodes
    elif wildcards.quantifier == 'starsolo':
        return rules.starsolo_quant.output.barcodes
    else:
        pass

def get_singlecell_bam(wildcards):
    if wildcards.quantifier == 'cellranger':
        return rules.cellranger_quant.output.bam
    elif wildcards.quantifier in ['starsolo', 'alevin']:
        return rules.starsolo_bam.output
    else:
        pass

def get_donor_vcf(wildcards):
    donor_dir = config['quant'].get('demultiplex', {}).get('donor_dir')
    if donor_dir is None:
        raise ValueError("missing configuration value: quant:demultiplex:donor_dir")
    if config['db']['reference_db'] == "ensembl":
        vcf_fn = "cellSNP.cells.vcf"
    else:
        vcf_fn = "cellSNP.cells.chr.vcf"
    return join(donor_dir, wildcards.sample, vcf_fn)
                
            
rule bam_rename_chromosomes:
    input:
        bam = get_singlecell_bam,
        mapping_file = '/mnt/archive/ext_cache/ext/ensembl/release-109/homo_sapiens/GRCh38/anno/ensembl_to_10xgenomics_2020-A_map.txt'
    output:
        bam = '{sample}_fixed.bam'
    container:
        'docker://lindenb/jvarkit:1b2aedf24'
    shell:
        'java -jar /opt/jvarkit/dist/jvarkit.jar bamrenamechr '
        '-f {input.mapping_file} '
        '--samoutputformat BAM '
        '--out {output.bam} '
        '{input.bam} '

rule sorted_vcf:
    input:
        vcf = join(REF_DIR, 'anno', 'common_variants.vcf'),
        bam = get_singlecell_bam
    output:
        vcf = join(DEMUX_DIR, 'common_variants_sorted.vcf.gz'),
        vcf2 = join(DEMUX_DIR, 'common_variants_sorted.vcf')
    params:
        script = src_gcf('scripts/sort_vcf_as_bam.py')
    threads:
        4
    container:
        'docker://' + config['docker']['vcftools']
    shell:
        'python {params.script} '
        '--vcf-file {input.vcf} '
        '--bam-file {input.bam} '
        '--output-file {output.vcf} '
        
rule freemuxlet_pileup:
    input:
        bam = get_singlecell_bam,
        barcodes = get_singlecell_barcodes,
        vcf = join(DEMUX_DIR, 'common_variants_sorted.vcf'),
    output:
        cel = join(DEMUX_DIR,  'freemuxlet', 'pileup.cel.gz')
    params:
        freemuxlet_dir = join(DEMUX_DIR,  'freemuxlet')
    threads:
        24
    container:
        'docker://' + config['docker']['popscle']
    shell:
        'popscle dsc-pileup '
        '--sam {input.bam} '
        '--vcf {input.vcf} '
        '--group-list {input.barcodes} '
        '--out {params.freemuxlet_dir}/pileup '

rule cellsnp_pileup_1a: # pileup with defined snps and barcodes (singlecell)
    input:
        bam = get_singlecell_bam,
        barcodes = get_singlecell_barcodes,
        vcf = join(REF_DIR, 'anno', 'common_variants.vcf')
    output:
        vcf = join(DEMUX_DIR,  'cellsnp', 'cellSNP.base.vcf.gz'),
        samples = join(DEMUX_DIR,  'cellsnp', 'cellSNP.samples.tsv'),
        mtx_ad = join(DEMUX_DIR,  'cellsnp', 'cellSNP.tag.AD.mtx'),
        mtx_dp = join(DEMUX_DIR,  'cellsnp', 'cellSNP.tag.DP.mtx'),
        mtx_other = join(DEMUX_DIR,  'cellsnp', 'cellSNP.tag.OTH.mtx')
    params:
        cellsnp_dir = join(DEMUX_DIR,  'cellsnp'),
        input_bam = lambda wildcards, input: input.bam if isinstance(input.bam, str) else ','.join(input.bam),
        minMAF = 0.1,
        minCOUNT = 20,
        args = '--gzip '
    threads:
        16
    container:
        'docker://' + config['docker']['cellsnp-lite']
    shell:
        'cellsnp-lite '
        '-s {params.input_bam} '
        '-b {input.barcodes} '
        '-R {input.vcf} '
        '-O {params.cellsnp_dir} '
        '--minMAF {params.minMAF} '
        '--minCOUNT {params.minCOUNT} ' 
        ' {params.args} '
        '--nproc {threads} '

rule vireo_noref:
    input:
        rules.cellsnp_pileup_1a.output
    output:
        donor_ids = join(DEMUX_DIR,  'vireo', 'donor_ids.tsv'),
        doublet = join(DEMUX_DIR,  'vireo', 'prob_doublet.tsv.gz'),
        singlet = join(DEMUX_DIR,  'vireo', 'prob_singlet.tsv.gz'),
        summary = join(DEMUX_DIR,  'vireo', 'summary.tsv')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        vireo_dir = join(DEMUX_DIR,  'vireo'),
        cellsnp_dir= rules.cellsnp_pileup_1a.params.cellsnp_dir
    container:
        'docker://' + config['docker']['vireo']
    shell:
        'vireo '
        '-c {params.cellsnp_dir} '
        '-o {params.vireo_dir} '
        '-N {params.n} '

rule vireo_ref:
    input:
        donor_vcf = get_donor_vcf,
        single_cell_vcf = rules.cellsnp_pileup_1a.output
    output:
        donor_ids = join(DEMUX_DIR, 'with_reference', 'vireo', 'donor_ids.tsv'),
        doublet = join(DEMUX_DIR, 'with_reference', 'vireo', 'prob_doublet.tsv.gz'),
        singlet = join(DEMUX_DIR, 'with_reference', 'vireo', 'prob_singlet.tsv.gz'),
        summary = join(DEMUX_DIR, 'with_reference', 'vireo', 'summary.tsv')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        vireo_dir = join(DEMUX_DIR,  'with_reference', 'vireo'),
        cellsnp_dir= rules.cellsnp_pileup_1a.params.cellsnp_dir
    threads:
        12
    container:
        'docker://' + config['docker']['vireo']
    shell:
        'vireo '
        '-d {input.donor_vcf} '
        '-c {params.cellsnp_dir} '
        '-o {params.vireo_dir} '
        '-N {params.n} '
        '--callAmbientRNAs '
        '--nproc {threads} '

rule vireo_droplet_type:
    input:
        join('{anything}', 'vireo', 'donor_ids.tsv')
    output:
        join('{anything}', 'vireo', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/vireo_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} -i {input} -o {output} '

rule souporcell_noref:
    input:
        bam = get_singlecell_bam,
        barcodes = get_singlecell_barcodes,
        vcf = join(REF_DIR, 'anno', 'common_variants.vcf'),
        genome = join(REF_DIR, 'fasta', 'genome.fa')
    output:
        vcf = join(DEMUX_DIR,  'souporcell', 'cluster_genotypes.vcf'),
        donor_ids = join(DEMUX_DIR, 'souporcell', 'clusters.tsv')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        souporcell_dir = join(DEMUX_DIR,  'souporcell'),
        args = '--ploidy 2 --min_alt 10 --min_ref 10 --max_loci 2048 --skip_remap True '
    threads:
        24
    container:
        'shub://' + config['docker']['souporcell']
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
        
rule souporcell_ref:
    input:
        bam = get_singlecell_bam,
        barcodes = get_singlecell_barcodes,
        donor_vcf = get_donor_vcf,
        genome = join(REF_DIR, 'fasta', 'genome.fa'),
    output:
        vcf = join(DEMUX_DIR,  'with_reference', 'souporcell', 'cluster_genotypes.vcf'),
        donor_ids = join(DEMUX_DIR,  'with_reference', 'souporcell', 'clusters.tsv')
    params:
        n = config['quant'].get('demultiplex', {}).get('n_individuals', 4),
        souporcell_dir = join(DEMUX_DIR,  'with_reference', 'souporcell'),
        #donor_ids = lambda wildcards: ' '.join(get_donor_ids(wildcards)),
        args = '--skip_remap True '
    threads:
        24
    container:
        'shub://' + config['docker']['souporcell']
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
        #'--known_genotypes_sample_names {params.donor_ids} '

rule souporcell_droplet_type:
    input:
        join('{anything}', 'souporcell', 'clusters.tsv')
    output:
        join('{anything}', 'souporcell', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/souporcell_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} -i {input} -o {output} '

rule freemuxlet_noref:
    input:
        barcodes = get_singlecell_barcodes,
        cel = join(DEMUX_DIR,  'freemuxlet', 'pileup.cel.gz')
    output:
        vcf = join(DEMUX_DIR,  'freemuxlet', 'freemuxlet.clust1.vcf.gz'),
        donor_ids = join(DEMUX_DIR,  'freemuxlet', 'freemuxlet.clust1.samples.gz')
    params:
        freemuxlet_dir = join(DEMUX_DIR,  'freemuxlet'),
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

rule freemuxlet_droplet_type:
    input:
        join('{anything}', 'freemuxlet', 'clusters.tsv')
    output:
        join('{anything}', 'freemuxlet', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/freemuxlet_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} -i {input} -o {output} '

rule demuxlet_ref:
    input:
        barcodes = get_singlecell_barcodes,
        cel = join(DEMUX_DIR,  'freemuxlet', 'pileup.cel.gz'),
        donor_vcf = get_donor_vcf
    output:
        donor_ids = join(DEMUX_DIR,  'with_reference', 'demuxlet', 'demuxlet.best')
    params:
        freemuxlet_dir = join(DEMUX_DIR,  'freemuxlet'),
        demuxlet_prefix = join(DEMUX_DIR,  'with_reference', 'demuxlet', 'demuxlet'),
        demuxlet_dir = join(DEMUX_DIR,  'with_reference', 'demuxlet'),
        field = 'GT',
        args = ' '
    threads:
        24
    container:
        'docker://' + config['docker']['popscle']
    shell:
        'popscle demuxlet '
        '--plp {params.freemuxlet_dir}/pileup '
        '--vcf {input.donor_vcf} '
        '--field {params.field} '
        '--group-list {input.barcodes} '
        '--out {params.demuxlet_prefix} '
        '{params.args} '

rule demuxlet_droplet_type:
    input:
        join('{anything}', 'demuxlet', 'demuxlet.best')
    output:
        join('{anything}', 'demuxlet', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/freemuxlet_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} -i {input} -o {output} '

rule demuxalot_noref:
    input:
        barcodes = get_singlecell_barcodes,
        bam = get_singlecell_bam,
        cluster_vcf = join(DEMUX_DIR,  'freemuxlet', 'freemuxlet.clust1.vcf.gz')
    output:
        join(DEMUX_DIR,  'demuxalot', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/run_demuxalot.py')
    threads:
        12
    container:
        'docker://' + config['docker']['demuxalot']
    shell:
        'python {params.script} '
        '-b {input.barcodes} '
        '-i {input.bam} '
        '-v {input.cluster_vcf} '
        '-o {output} '

rule demuxalot_ref:
    input:
        barcodes = get_singlecell_barcodes,
        bam = get_singlecell_bam,
        donor_vcf = get_donor_vcf
    output:
        join(DEMUX_DIR,  'with_reference', 'demuxalot', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/run_demuxalot.py')
    threads:
        12
    container:
        'docker://' + config['docker']['demuxalot']
    shell:
        'python {params.script} '
        '-b {input.barcodes} '
        '-i {input.bam} '
        '-v {input.donor_vcf} '
        '--threads {threads} '
        '-o {output} '


def _get_demuxafy_demux_methods():
    """returns the best combo of doublet detection methods according to demuxafy when multiplexing is available."""
    if not SAMPLE_MULTIPLEXING:
        raise ValueError
    n = config['quant']['demultiplex'].get('n_reference_samples', 4)
    if config['quant']['demultiplex']['with_reference']:
        if n < 5:
            doublet_methods = ['scds', 'scdblfinder', 'scrublet']
            demux_methods = ['with_reference/demuxalot', 'with_reference/vireo']
        elif n <= 20:
            doublet_methods = ['scds', 'scdblfinder', 'doublet_detection']
            demux_methods = ['with_reference/demuxalot', 'with_reference/dropulation']
        else:
            doublet_methods = []
            demux_methods = ['with_reference/dropulation', 'with_reference/demuxalot']
    else: # no reference
        if n < 5:
            doublet_method = ['scds', 'scdblfinder', 'scrublet']
            demux_methods = ['soupercell', 'vireo']
        elif n <= 20:
            doublet_methods = ['scds', 'scdblfinder', 'doublet_detection']
            demux_methods = ['soupercell', 'vireo']
        else:
            doublet_methods = []
            demux_methods = ['vireo']

    return doublet_methods, demux_methods

def _get_default_demux_methods(mtype='doublet'):
    """overrides method in doublets.smk
    """
    if config['quant']['demultiplex']['with_reference']:
        doublet_methods = ['socube', 'scdblfinder']
        demux_methods = ['with_reference/demuxalot',]
    else:
        doublet_methods = ['socube', 'scdblfinder']
        demux_methods = ['vireo']
    if mtype == 'doublet':
        return doublet_methods
    elif mtype == 'demux':
        return demux_methods
    
def get_multiplex_methods(test_all=False):
    if test_all:
        # reference souporcell failing (maybe donor-id related?, works on no-ref)
        doublet_methods = ['scds', 'solo', 'scrublet', 'doubletdetection', 'scdblfinder', 'socube']
        demux_methods = ['with_reference/vireo', 'with_reference/demuxalot', 'with_reference/demuxlet']
    else:
        demux_methods = config['quant'].get('demultiplex', {}).get('method')
        doublet_methods = config['quant'].get('doublet_detection', {}).get('method')
        
        if demux_methods is None or demux_methods == 'default':
            demux_methods =  _get_default_demux_methods('demux')
        elif demux_methods == 'demuxafy':
            doublet_methods, demux_methods = _get_demuxafy_demux_methods()
        elif demux_methods[0] == 'skip':
            raise NotImplementedError        
        else:
            demux_methods = demux_methods.split(',')
            if config['quant']['demultiplex']['with_reference']:
                demux_methods = [os.path.join("with_reference", i) for i in demux_methods]
            
        if doublet_methods is None or doublet_methods == 'default':
            doublet_methods =  _get_default_demux_methods('doublet')
        else:
            doublet_methods = doublet_methods.split(',')
        

    dbl_files = expand(join(QUANT_INTERIM, '{{quantifier}}', '{{sample}}', 'doublets', '{method}', 'doublet_type.tsv'),
                       method=doublet_methods)
    demux_files = expand(join(QUANT_INTERIM, '{{quantifier}}', '{{sample}}', 'demultiplexing', '{method}', 'droplet_type.tsv'),
                         method=demux_methods)
    return [dbl_files, demux_files]
