# bfq_level2_common.smk
# Shared helpers for all level2 collectors

BFQ_LEVEL2_ALL = []

# Expect these symbols from the parent workflow:
# - QUANT_INTERIM, BFQ_INTERIM, CR_INTERIM
# - AGGR_IDS, SAMPLES
# - config, CB_OUTPUT
# - BFQ_ALL (final fan-in list)

from os.path import join
from snakemake.shell import shell

def symlink(src, dst):
    shell("ln -sfnr {src} {dst}")

def notebook_inputs(method):
    return [expand(join(QUANT_INTERIM, "aggregate", method, "scanpy", "notebooks", "{aggr_id}_pp.html"), aggr_id=AGGR_IDS),
            expand(join(QUANT_INTERIM, "aggregate", method, "scanpy", "notebooks", "{aggr_id}_pp.ipynb"), aggr_id=AGGR_IDS)]

def notebook_outputs():
    return [expand(join(BFQ_INTERIM, "notebooks", "{aggr_id}_preprocess.html"), aggr_id=AGGR_IDS),
            expand(join(BFQ_INTERIM, "notebooks", "{aggr_id}_preprocess.ipynb"), aggr_id=AGGR_IDS)]

def exprs_aggr_input(method):
    return expand(join(QUANT_INTERIM, "aggregate", method, "scanpy", "{aggr_id}_filtered.h5ad"), aggr_id=AGGR_IDS)

def exprs_aggr_output():
    return expand(join(BFQ_INTERIM, "exprs", "scanpy", "{aggr_id}_filtered.h5ad"), aggr_id=AGGR_IDS)


rule bfq_level2_umap_png:
    input:
        join(BFQ_INTERIM, "exprs", "scanpy", "all_samples_filtered.h5ad")
    output:
        join(BFQ_INTERIM, 'figs', 'umap_all_samples_mqc.png')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -o {output}'
            
rule bfq_level2_umap_yaml:
    input:
        join(BFQ_INTERIM, "exprs", "scanpy", "all_samples_filtered.h5ad")
    output:
        join(BFQ_INTERIM, 'figs', 'all_samples_mqc.yaml')
    params:
        script = src_gcf('scripts/plotpca.py')
    container:
        'docker://' + config['docker']['scanpy']
    shell:
        'python {params.script} {input} -o {output}'
