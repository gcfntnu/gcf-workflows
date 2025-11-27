#-*- mode:snakemake -*-
from snakemake.utils import min_version
min_version("6.0")


include:
    "bfq_level2_common.smk"

m = config['quant']['method']
if m == "cellranger":
    include: "bfq_level2_cellranger.smk"

elif m == "splitpipe":
    include:
        "bfq_level2_splitpipe.smk"

elif m == "parsebio_starsolo":
    include:
        "bfq_level2_parsebio_starsolo.smk"

else:
    pass
    
BFQ_ALL.extend(BFQ_LEVEL2_ALL)

