rule multiqc_config:
    input:
        header_template = srcdir(join('misc', 'multiqc', 'mqc_header.txt')),
        config_template = srcdir(join('misc', 'multiqc', 'multiqc_config-{}.yaml'.format(WORKFLOW))),
        sample_info = join(INTERIM_DIR, 'sample_info.tsv'),
    output:
        mqc_config = join(BFQ_INTERIM, '.multiqc_config.yaml')
    params:
        script = srcdir('misc/multiqc/create_mqc_config.py'),
        org = config['organism'],
        project_id = config['project_id'],
        machine = config['machine'],
        read_geometry = ','.join([str(x) for x in config['read_geometry']]),
        libprep = config['libprepkit'],
        repo_dir = srcdir(os.path.dirname('main.config')),
    singularity:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '-p {params.project_id} '
        '-S {input.sample_info} '
        '--organism {params.org} '
        '--libkit "{params.libprep}" '
        '--machine "{params.machine}" '
        '--read-geometry {params.read_geometry} '
        '--repo-dir {params.repo_dir} '
        '--header-template {input.header_template} '
        '--config-template {input.config_template} '
        '-o {output.mqc_config} '

def get_mqc_modules():
    if WORKFLOW in config['multiqc'].keys():
        modules = config['multiqc'][WORKFLOW].get('modules', "").split(',')
        return ('-m ' + ' -m '.join(modules)) if modules else ""
    else:
        return ""

rule multiqc_report:
    input:
        bfq = BFQ_ALL,
        mqc_config = rules.multiqc_config.output.mqc_config,
    output:
        report = join(BFQ_INTERIM, 'multiqc_{}.html'.format('_'.join(config['project_id']))),
    params:
        modules = get_mqc_modules(),
        extra_args = '-f -q --interactive ',
        analysis_directory = BFQ_INTERIM,
    singularity:
        'docker://' + config['docker']['multiqc']
    shell:
        'multiqc '
        '--config {input.mqc_config} '
        '{params.modules} '
        '--filename {output.report} '
        '{params.extra_args} '
        '{params.analysis_directory} '


