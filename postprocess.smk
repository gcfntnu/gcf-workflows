
PROJECT_ID = config['project_id']
PROJECT_ID = PROJECT_ID[0] if isinstance(PROJECT_ID, (list, tuple)) else PROJECT_ID


# report on read geometry before any filter steps
_original_read_geometry = config['read_geometry']
if 'delta_readlen' in config:
    for i, val in enumerate(config['delta_readlen']):
        _original_read_geometry[i] = int(config['read_geometry'][i]) - int(val)

rule download_geo_template:
    output:
        temp('seq_template.xlsx')
    params:
        proxy = config.get('proxy', {}).get('wget', ''),
        url = 'https://www.ncbi.nlm.nih.gov/geo/info/examples/seq_template.xlsx'
    shell:
        'wget {params.proxy} {params.url} -O- > {output} '

rule geo_md5sum_sample:
    input:
        unpack(get_merged_fastq)
    params:
        fastq_dir = lambda wildcards, input: os.path.dirname(input.R1)
    output:
        temp('{sample}.md5')
    shell:
        'md5sum {params.fastq_dir}/{wildcards.sample}*fastq.gz > {output}'
        
rule geo_md5sum:
    input:
        expand(rules.geo_md5sum_sample.output, sample=SAMPLES)
    output:
        temp('fastq.md5')
    shell:
        'cat {input} > {output}'

        
rule geo_template:
    input:
        bfq = BFQ_ALL,
        excel_template =  rules.download_geo_template.output,
        sample_info = join(INTERIM_DIR, 'sample_info.tsv'),
        processed_data = GEO_PROCESSED_FILES,
        fastq_md5 = rules.geo_md5sum.output
    output:
        geo_filled_template = join(BFQ_INTERIM, 'data_submission', 'geo', 'seq_template_' + PROJECT_ID + '.xlsx')
    params:
        script = srcdir('misc/protocols/geoguessr.py'),
        read_geometry = ','.join([str(x) for x in _original_read_geometry]),
        pep = workflow.pepfile,
        assembly = config['db'].get(config['db'].get('reference_db'), {}).get('assembly', 'NA'),
        repo_dir = srcdir(os.path.dirname('main.config')),
    singularity:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '-S {input.sample_info} '
        '--template {input.excel_template} '
        '--md5sum-file {input.fastq_md5} '
        '--processed-data {input.processed_data} '
        '--read-geometry {params.read_geometry} '
        '--pep {params.pep} '
        '--assembly {params.assembly} '
        '--repo-dir {params.repo_dir}  '
        '-o {output.geo_filled_template} '

rule multiqc_config:
    input:
        header_template = srcdir(join('misc', 'multiqc', 'mqc_header.txt')),
        config_template = srcdir(join('misc', 'multiqc', 'multiqc_config-{}.yaml'.format(WORKFLOW))),
        sample_info = join(INTERIM_DIR, 'sample_info.tsv'),
    output:
        mqc_config = join(BFQ_INTERIM, '.multiqc_config.yaml')
    params:
        script = srcdir('misc/multiqc/create_mqc_config.py'),
        project_id = PROJECT_ID,
        machine = config['machine'],
        read_geometry = ','.join([str(x) for x in _original_read_geometry]),
        libprep = config['libprepkit'],
        repo_dir = srcdir(os.path.dirname('main.config')),
        pep = workflow.pepfile
    singularity:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '-p {params.project_id} '
        '-S {input.sample_info} '
        '--libkit "{params.libprep}" '
        '--machine "{params.machine}" '
        '--read-geometry {params.read_geometry} '
        '--repo-dir {params.repo_dir} '
        '--header-template {input.header_template} '
        '--config-template {input.config_template} '
        '--pep {params.pep} '
        '-o {output.mqc_config} '

def get_mqc_modules():
    if WORKFLOW in config['multiqc'].keys():
        modules = config['multiqc'][WORKFLOW].get('modules', "").split(',')
        return ('-m ' + ' -m '.join(modules)) if modules else ""
    else:
        return ""

def bfq_search_paths(*args):
    """run time check for multiqc search paths in bfq directory
    """
    valid_paths = []
    log_path = join(BFQ_INTERIM, 'logs')
    fig_path = join(BFQ_INTERIM, 'figs')
    if os.path.exists(log_path):
        valid_paths.append(log_path)
    if os.path.exists(fig_path):
        valid_paths.append(fig_path)
    if len(valid_paths) == 0:
        valid_paths = BFQ_INTERIM
    return valid_paths
    
rule multiqc_report:
    input:
        bfq = BFQ_ALL,
        mqc_config = rules.multiqc_config.output.mqc_config,
    output:
        report = join(BFQ_INTERIM, 'multiqc_{}.html'.format(PROJECT_ID)),
    params:
        modules = get_mqc_modules(),
        extra_args = '-f -q --interactive ',
        search_paths = lambda wildcards: bfq_search_paths()
    singularity:
        'docker://' + config['docker']['multiqc']
    shell:
        'multiqc '
        '--config {input.mqc_config} '
        '{params.modules} '
        '--filename {output.report} '
        '{params.extra_args} '
        '{params.search_paths} '


