"""demultiplexing by natural variation

note: donor vcf files are assumed to contain the specifc samples to be demultiplexed and no more. 
"""
DEMUX_DIR = join(QUANT_INTERIM, '{quantifier}', '{sample}', 'demultiplexing')
DEMUX_CONFIG = config['quant']['demultiplex']

CELLSNP_CONFIG = DEMUX_CONFIG['cellsnp']
DEMUXALOT_CONFIG = DEMUX_CONFIG['demuxalot']
VIREO_CONFIG = DEMUX_CONFIG['vireo']
SOUPORCELL_CONFIG = DEMUX_CONFIG['souporcell']

DEMUX_METHODS = [m.strip() for m in DEMUX_CONFIG['method'].split(',')]


def demux_method_uses_reference(method):
    return method.endswith('_ref') or method in ['vireo_ref_rescue', 'souporcell_ref_demuxafy']


if DEMUX_CONFIG['method'] != 'skip':
    if DEMUX_CONFIG['n_individuals'] < 1:
        raise ValueError(
            'quant.demultiplex.n_individuals must be set when '
            'demultiplexing is enabled.'
        )

    if any(demux_method_uses_reference(m) for m in DEMUX_METHODS) and not DEMUX_CONFIG.get('donor_dir'):
        raise ValueError(
            'quant.demultiplex.donor_dir must be set when '
            'reference-based demultiplexing is enabled.'
        )


def get_singlecell_barcodes(wildcards):
    if wildcards.quantifier == 'cellranger':
        if config['quant'].get('cellbender_filter', False):
            return rules.cellranger_cellbender.output.aggr
        return rules.cellranger_quant.output.filt_barcodes
    elif wildcards.quantifier == '10x_starsolo':
        return rules.starsolo_quant.output.barcodes
    else:
        raise ValueError

def get_singlecell_bam(wildcards):
    if wildcards.quantifier == 'cellranger':
        return rules.cellranger_quant.output.bam
    elif wildcards.quantifier in ['10x_starsolo', 'alevin']:
        return rules.starsolo_bam.output
    else:
        raise ValueError

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
        cel = join(DEMUX_DIR,  'freemuxlet_noref', 'pileup.cel.gz')
    params:
        freemuxlet_dir = join(DEMUX_DIR,  'freemuxlet_noref')
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
        min_maf = CELLSNP_CONFIG["min_maf"],
        min_count = CELLSNP_CONFIG["min_count"],
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
        '--minMAF {params.min_maf} '
        '--minCOUNT {params.min_count} ' 
        ' {params.args} '
        '--nproc {threads} '

rule vireo_donor_subset:
    input:
        donor_vcf = get_donor_vcf,
        cellsnp_vcf = rules.cellsnp_pileup_1a.output.vcf,
    output:
        vcf = join(DEMUX_DIR, "vireo_ref", "donor_subset.vcf"),
    container:
        "docker://biocontainers/bcftools:v1.9-1-deb_cv1"
    shell:
        'bcftools view '
        '{input.donor_vcf} '
        '-T {input.cellsnp_vcf} '
        '-Ov '
        '-o {output.vcf}'

rule vireo_noref:
    input:
        rules.cellsnp_pileup_1a.output
    output:
        donor_ids = join(DEMUX_DIR, 'vireo_noref', 'donor_ids.tsv'),
        doublet = join(DEMUX_DIR, 'vireo_noref', 'prob_doublet.tsv.gz'),
        singlet = join(DEMUX_DIR, 'vireo_noref', 'prob_singlet.tsv.gz'),
        summary = join(DEMUX_DIR, 'vireo_noref', 'summary.tsv')
    params:
        n = DEMUX_CONFIG["n_individuals"],
        vireo_dir = join(DEMUX_DIR, 'vireo_noref'),
        cellsnp_dir = rules.cellsnp_pileup_1a.params.cellsnp_dir,
        rand_seed = VIREO_CONFIG.get("rand_seed", 1)
    threads:
        24
    container:
        'docker://' + config['docker']['vireo']
    shell:
        'vireo '
        '-c {params.cellsnp_dir} '
        '-o {params.vireo_dir} '
        '-N {params.n} '
        '--nproc {threads} '
        '--randSeed {params.rand_seed}'

rule vireo_ref:
    input:
        donor_vcf = rules.vireo_donor_subset.output.vcf,
        single_cell_vcf = rules.cellsnp_pileup_1a.output,
    output:
        donor_ids = join(DEMUX_DIR, "vireo_ref", "donor_ids.tsv"),
        doublet = join(DEMUX_DIR, "vireo_ref", "prob_doublet.tsv.gz"),
        singlet = join(DEMUX_DIR, "vireo_ref", "prob_singlet.tsv.gz"),
        summary = join(DEMUX_DIR, "vireo_ref", "summary.tsv"),
    params:
        n = DEMUX_CONFIG["n_individuals"],
        vireo_dir = join(DEMUX_DIR, "vireo_ref"),
        cellsnp_dir = rules.cellsnp_pileup_1a.params.cellsnp_dir,
        force_learn_gt = "--forceLearnGT" if VIREO_CONFIG.get("force_learn_gt", False) else ""
    threads:
        24
    container:
        "docker://" + config["docker"]["vireo"]
    shell:
        'vireo '
        '-d {input.donor_vcf} '
        '-c {params.cellsnp_dir} '
        '-o {params.vireo_dir} '
        '-N {params.n} '
        '--nproc {threads} '
        '{params.force_learn_gt}'


rule vireo_droplet_type:
    input:
        donor_ids = join('{anything}', 'vireo_{mode}', 'donor_ids.tsv')
    output:
        droplet_type = join('{anything}', 'vireo_{mode}', 'droplet_type.tsv')
    wildcard_constraints:
        mode = 'noref|ref'
    params:
        script = src_gcf('scripts/vireo_summary.py'),
        min_prob = VIREO_CONFIG["min_prob"],
        min_doublet_prob = VIREO_CONFIG["min_doublet_prob"],
        min_vars = VIREO_CONFIG["min_vars"],
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '-i {input.donor_ids} '
        '-o {output.droplet_type} '
        '--min-prob {params.min_prob} '
        '--min-doublet-prob {params.min_doublet_prob} '
        '--min-vars {params.min_vars} '


rule souporcell_noref:
    input:
        bam = get_singlecell_bam,
        barcodes = get_singlecell_barcodes,
        vcf = join(REF_DIR, 'anno', 'common_variants.vcf'),
        genome = join(REF_DIR, 'fasta', 'genome.fa')
    output:
        vcf = join(DEMUX_DIR,  'souporcell_noref', 'cluster_genotypes.vcf'),
        donor_ids = join(DEMUX_DIR, 'souporcell_noref', 'clusters.tsv')
    params:
        n = DEMUX_CONFIG["n_individuals"],
        souporcell_dir = join(DEMUX_DIR,  'souporcell_noref'),
        skip_remap = '--skip_remap True' if SOUPORCELL_CONFIG['skip_remap'] else ''
    threads:
        24
    container:
        'shub://' + config['docker']['souporcell']
    shell:
        'mkdir -p {params.souporcell_dir}/logs && '
        'souporcell_pipeline.py '
        '-i {input.bam} '
        '-b {input.barcodes} '
        '-f {input.genome} '
        '--common_variants {input.vcf} '
        '-t {threads} '
        '-o {params.souporcell_dir} '
        '-k {params.n} '
        '{params.skip_remap} '
        

rule souporcell_ref:
    input:
        bam = get_singlecell_bam,
        barcodes = get_singlecell_barcodes,
        donor_vcf = get_donor_vcf,
        genome = join(REF_DIR, "fasta", "genome.fa")
    output:
        vcf = join(DEMUX_DIR, "souporcell_ref", "cluster_genotypes.vcf"),
        donor_ids = join(DEMUX_DIR, "souporcell_ref", "clusters.tsv")
    params:
        n = DEMUX_CONFIG["n_individuals"],
        souporcell_dir = join(DEMUX_DIR, "souporcell_ref"),
        skip_remap = '--skip_remap True' if SOUPORCELL_CONFIG['skip_remap'] else ''
    threads:
        24
    container:
        "library://wheaton5/souporcell/souporcell:release"
    shell:
        'mkdir -p {params.souporcell_dir}/logs && '
        'souporcell_pipeline.py '
        '-i {input.bam} '
        '-b {input.barcodes} '
        '-f {input.genome} '
        '--known_genotypes {input.donor_vcf} '
        '-t {threads} '
        '-o {params.souporcell_dir} '
        '-k {params.n} '
        '{params.skip_remap} '


rule souporcell_ref_demuxafy:
    input:
        donor_vcf = get_donor_vcf,
        cluster_vcf = rules.souporcell_noref.output.vcf,
        clusters = rules.souporcell_noref.output.donor_ids
    output:
        genotype_key = join(DEMUX_DIR, 'souporcell_ref_demuxafy', 'Genotype_ID_key.txt'),
        correlations = join(DEMUX_DIR, 'souporcell_ref_demuxafy', 'ref_clust_pearson_correlations.tsv'),
        droplet_type = join(DEMUX_DIR, 'souporcell_ref_demuxafy', 'droplet_type.tsv')
    params:
        outdir = join(DEMUX_DIR, 'souporcell_ref_demuxafy')
    container:
        'docker://' + config['docker']['demuxafy-match']
    shell:
        'Assign_Indiv_by_Geno.R '
        '-r {input.donor_vcf} '
        '-c {input.cluster_vcf} '
        '-o {params.outdir} '
        '--clusters {input.clusters} '


rule souporcell_droplet_type:
    input:
        join('{anything}', 'souporcell_{mode}', 'clusters.tsv')
    output:
        join('{anything}', 'souporcell_{mode}', 'droplet_type.tsv')
    wildcard_constraints:
        mode = 'noref|ref'
    params:
        script = src_gcf('scripts/souporcell_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} -i {input} -o {output} '

rule freemuxlet_noref:
    input:
        barcodes = get_singlecell_barcodes,
        cel = join(DEMUX_DIR,  'freemuxlet_noref', 'pileup.cel.gz')
    output:
        vcf = join(DEMUX_DIR,  'freemuxlet_noref', 'freemuxlet.clust1.vcf.gz'),
        donor_ids = join(DEMUX_DIR,  'freemuxlet_noref', 'freemuxlet.clust1.samples.gz')
    params:
        freemuxlet_dir = join(DEMUX_DIR,  'freemuxlet_noref'),
        n = DEMUX_CONFIG["n_individuals"]
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
        join('{anything}', 'freemuxlet_noref', 'freemuxlet.clust1.samples.gz')
    output:
        join('{anything}', 'freemuxlet_noref', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/freemuxlet_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} -i {input} -o {output} '

rule demuxlet_ref:
    input:
        barcodes = get_singlecell_barcodes,
        cel = join(DEMUX_DIR,  'freemuxlet_noref', 'pileup.cel.gz'),
        donor_vcf = get_donor_vcf
    output:
        donor_ids = join(DEMUX_DIR,  'demuxlet_ref', 'demuxlet.best')
    params:
        freemuxlet_dir = join(DEMUX_DIR,  'freemuxlet_noref'),
        demuxlet_prefix = join(DEMUX_DIR,  'demuxlet_ref', 'demuxlet'),
        demuxlet_dir = join(DEMUX_DIR,  'demuxlet_ref'),
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
        join('{anything}', 'demuxlet_ref', 'demuxlet.best')
    output:
        join('{anything}', 'demuxlet_ref', 'droplet_type.tsv')
    params:
        script = src_gcf('scripts/demuxlet_summary.py')
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} -i {input} -o {output} '

rule demuxalot_noref:
    input:
        barcodes = get_singlecell_barcodes,
        bam = get_singlecell_bam,
        cluster_vcf = join(DEMUX_DIR,  'freemuxlet_noref', 'freemuxlet.clust1.vcf.gz')
    output:
        droplet_type = join(DEMUX_DIR,  'demuxalot_noref', 'droplet_type.tsv'),
        snp_summary = join(DEMUX_DIR, 'demuxalot_noref', 'snp_summary.tsv'),
        snp_plot = join(DEMUX_DIR, 'demuxalot_noref', 'snp_summary.pdf')
    params:
        script = src_gcf('scripts/run_demuxalot.py'),
        genotype_mode = DEMUXALOT_CONFIG["genotype_mode"],
        doublet_prior = DEMUXALOT_CONFIG["doublet_prior"],
        min_prob = DEMUXALOT_CONFIG["min_prob"],
        vcf_prior_strength = DEMUXALOT_CONFIG["vcf_prior_strength"],
    threads:
        12
    container:
        'docker://' + config['docker']['demuxalot']
    shell:
        'python {params.script} '
        '-b {input.barcodes} '
        '-i {input.bam} '
        '-v {input.cluster_vcf} '
        '--genotype-mode {params.genotype_mode} '
        '--doublet-prior {params.doublet_prior} '
        '--min-prob {params.min_prob} '
        '--vcf-prior-strength {params.vcf_prior_strength} '
        '--threads {threads} '
        '-o {output.droplet_type} '

rule demuxalot_ref:
    input:
        barcodes = get_singlecell_barcodes,
        bam = get_singlecell_bam,
        donor_vcf = get_donor_vcf
    output:
        droplet_type = join(DEMUX_DIR,  'demuxalot_ref', 'droplet_type.tsv'),
        snp_summary = join(DEMUX_DIR, 'demuxalot_ref', 'snp_summary.tsv'),
        snp_plot = join(DEMUX_DIR, 'demuxalot_ref', 'snp_summary.pdf'),
        learned_genotypes = join(DEMUX_DIR, 'demuxalot_ref', 'genotypes.pq')
    params:
        script = src_gcf('scripts/run_demuxalot.py'),
        genotype_mode = DEMUXALOT_CONFIG["genotype_mode"],
        doublet_prior = DEMUXALOT_CONFIG["doublet_prior"],
        min_prob = DEMUXALOT_CONFIG["min_prob"],
        vcf_prior_strength = DEMUXALOT_CONFIG["vcf_prior_strength"],
    threads:
        12
    container:
        'docker://' + config['docker']['demuxalot']
        
    shell:
        'python {params.script} '
        '-b {input.barcodes} '
        '-i {input.bam} '
        '-v {input.donor_vcf} '
        '--genotype-mode {params.genotype_mode} '
        '--doublet-prior {params.doublet_prior} '
        '--min-prob {params.min_prob} '
        '--vcf-prior-strength {params.vcf_prior_strength} '
        '--learned-genotypes-output {output.learned_genotypes} '
        '--threads {threads} '
        '-o {output.droplet_type} '


def _get_demuxafy_demux_methods():
    """returns the best combo of doublet detection methods according to demuxafy when multiplexing is available."""
    if not SAMPLE_MULTIPLEXING:
        raise ValueError
    n = DEMUX_CONFIG["n_individuals"]
    if DEMUX_CONFIG.get('donor_dir'):
        if n < 5:
            doublet_methods = ['scds', 'scdblfinder', 'scrublet']
            demux_methods = ['vireo_ref']
        elif n <= 20:
            doublet_methods = ['scds', 'scdblfinder', 'doubletdetection']
            demux_methods = ['demuxalot_ref', 'vireo_ref']
        else:
            doublet_methods = []
            demux_methods = ['vireo_ref', 'demuxalot_ref']
    else: # no reference
        if n < 5:
            doublet_methods = ['scds', 'scdblfinder', 'scrublet']
            demux_methods = ['souporcell_noref', 'vireo_noref']
        elif n <= 20:
            doublet_methods = ['scds', 'scdblfinder', 'doubletdetection']
            demux_methods = ['souporcell_noref', 'vireo_noref']
        else:
            doublet_methods = []
            demux_methods = ['vireo_noref']

    return doublet_methods, demux_methods

def _get_default_demux_methods(mtype='doublet'):
    if DEMUX_CONFIG.get('donor_dir'):
        demux_methods = ['vireo_ref']
    else:
        demux_methods = ['vireo_noref']

    if mtype == 'doublet':
        return _get_default_methods()
    elif mtype == 'demux':
        return demux_methods

def get_multiplex_demux_methods():
    demux_methods = DEMUX_CONFIG.get('method')

    if demux_methods is None or demux_methods == 'default':
        return _get_default_demux_methods('demux')

    if demux_methods == 'demuxafy':
        _, demux_methods = _get_demuxafy_demux_methods()
        return demux_methods

    if demux_methods == 'skip':
        return []

    return [m.strip() for m in demux_methods.split(',')]


def get_multiplex_methods(test_all=False, n_cells=None):
    if test_all:
        doublet_methods = ['scds', 'solo', 'scrublet', 'doubletdetection', 'scdblfinder', 'socube']
        demux_methods = ['vireo_ref', 'demuxalot_ref', 'demuxlet_ref']
    else:
        demux_methods = get_multiplex_demux_methods()
        doublet_methods = config['quant'].get('doublet_detection', {}).get('method')

        if isinstance(doublet_methods, str):
            if doublet_methods == 'default':
                doublet_methods = _get_default_methods()
            elif doublet_methods == 'demuxafy':
                n_cells = n_cells or config['quant'].get('doublet_detection', {}).get('n_expected_cells')
                if n_cells is None:
                    raise ValueError("demuxafy requires n_cells argument")
                doublet_methods = _get_demuxafy_methods(n_cells=n_cells)
            else:
                doublet_methods = [m.strip() for m in doublet_methods.split(',')]
        elif doublet_methods is None:
            doublet_methods = _get_default_methods()

    dbl_files = expand(
        join(QUANT_INTERIM, '{{quantifier}}', '{{sample}}', 'doublets', '{method}', 'doublet_type.tsv'),
        method=doublet_methods
    )
    demux_files = expand(
        join(QUANT_INTERIM, '{{quantifier}}', '{{sample}}', 'demultiplexing', '{method}', 'droplet_type.tsv'),
        method=demux_methods
    )

    return [dbl_files, demux_files]

def multiplex_aggr_input(wildcards):
    inputs = {
        'droplet_types': [
            join(QUANT_INTERIM, wildcards.method, sample, 'demultiplexing', wildcards.multiplex_method, 'droplet_type.tsv')
            for sample in AGGR_IDS[wildcards.aggr_id]
        ]
    }

    if wildcards.method == 'cellranger' and AGGR_METHOD == 'cellranger':
        inputs['aggr_csv'] = join(QUANT_INTERIM, 'aggregate', 'description', f'{wildcards.aggr_id}_aggr.csv')

    return inputs


rule multiplex_droplet_type_aggr:
    input:
        unpack(multiplex_aggr_input)
    output:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'multiplexing', '{multiplex_method}', '{aggr_id}_droplet_type.tsv')
    params:
        script = src_gcf('scripts/aggr_barcode_info.py'),
        args = barcode_aggr_args
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '{input.droplet_types} '
        '{params.args} '
        '--output {output}'
