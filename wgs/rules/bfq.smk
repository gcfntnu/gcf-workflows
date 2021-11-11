
BFQ_ALL = []
include: "bfq/bfq_level1.smk"
include: "bfq/bfq_level2.smk"

rule bfq_all:
    input:
        BFQ_ALL
