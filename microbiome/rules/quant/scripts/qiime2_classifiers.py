"""Kit primers and exported classifier paths shared by the workflow and QIIME2."""

from os.path import join


def classifier_primers(libprep_config, libprep, taxonomy_db):
    primers_by_kit = {entry['libprep_name']: entry['db']['primers'] for entry in libprep_config.values() if 'primers' in entry.get('db', {})}
    if libprep not in primers_by_kit:
        raise ValueError('asked for: {}. Available libprep options: {}'.format(libprep, ', '.join(primers_by_kit)))
    if taxonomy_db not in ('silva', 'greengenes', 'unite'):
        raise ValueError('provided taxonomy-db {} is not supported!'.format(taxonomy_db))
    prefix = 'ITS' if taxonomy_db == 'unite' else 'V'
    primers = {region: pair for region, pair in primers_by_kit[libprep].items() if region.startswith(prefix)}
    if not primers:
        raise ValueError('Failed to identify correctly named regions (V* / ITS*)')
    return primers


def classifier_export_dir(classifier_dir, level, primer_pair):
    return join(classifier_dir, '{}_{}'.format(level, primer_pair))
