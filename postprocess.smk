
PROJECT_ID = config['project_id']
PROJECT_ID = PROJECT_ID[0] if isinstance(PROJECT_ID, (list, tuple)) else PROJECT_ID

# report on read geometry before any filter steps
_original_read_geometry = config['read_geometry']
if 'delta_readlen' in config:
    for i, val in enumerate(config['delta_readlen']):
        _original_read_geometry[i] = int(config['read_geometry'][i]) - int(val)
    
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
        pep = 'config.yaml' if config.get('skip_peppy', False) else workflow.pepfile
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


