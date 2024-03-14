#-*- mode:snakemake -*-

if not config['quant']['aggregate']['skip']:
    AGGR_IDS = collections.defaultdict(list)
    groupby = config['quant']['aggregate'].get('groupby', 'all_samples')
    for k, v in config['samples'].items():
        if groupby in v:
            aggr_id = v[groupby]
            AGGR_IDS[aggr_id].append(k)
        elif groupby == 'all_samples':
            AGGR_IDS['all_samples'].append(k)
        else:
            raise ValueError
        
include:
    'quant/cellranger.smk'
include:
    'quant/parse.smk'
include:
    'quant/alevin.smk'
include:
    'quant/umitools.smk'
include:
    'quant/star.smk'
include:
    'quant/velocyto.smk'
include:
    'quant/doublets.smk'
include:
    'qc/bam.smk'
    
QRULES = {'cellranger': rules.cellranger_quant.output,
          'starsolo': rules.starsolo_quant.output,
          'umitools': rules.umitools_quant.output,
          'alevin': rules.alevin_quant.output}
        
def get_quant(wildcards):
    method = config['quant'].get('method', 'cellranger')
    quant_rule = QRULES[method]
    if config['quant']['aggregate']['skip']:
        files = expand(quant_rule, samle=SAMPLES)
    else:
        aggr_method = config['quant']['aggregate'].get('method', 'cellranger')
        assert(aggr_method in ['cellranger', 'scanpy'])
        if aggr_method == 'cellranger':
            assert(method == aggr_method)
            files = expand(rules.cellranger_aggr.output, aggr_id=AGGR_IDS[wildcards.aggr_id])
        else:
            files = expand(rules.cellranger_aggr_scanpy.output, aggr_id=AGGR_IDS[wildcards.aggr_id])
    return files


def get_feature_info():
    """
    """
    return join(REF_DIR, 'anno', 'genes.tsv')

def get_barcode_info():
    if config['quant'].get('doublet_detection', {}).get('method') not in [None, 'skip']:
        #config['quant'].get('demultiplex', {}).get('method') not in [None, 'skip']
        doublets = join(QUANT_INTERIM, 'aggregate', config['quant']['method'] , 'all_samples_droplet_type.tsv')
        return doublets
    return None
