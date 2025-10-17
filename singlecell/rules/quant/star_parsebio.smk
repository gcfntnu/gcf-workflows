
#-*- mode:snakemake -*-
"""

"""
groupname = "_star_"

import tempfile
from typing import Iterable, Optional, Tuple, List


# Parse adapters
PRE_TSO_SEQ   = config["quant"].get("tso", "AACGCAGAGTGAATGGG")
LINKER_RC     = config["quant"].get("l21_rc", "AACGCAGAGTGAATGGG")
# Defaults from config 
STARSOLO_FEATURE = config["quant"].get("starsolo", {}).get("feature_count", "GeneFull_Ex50pAS")
STARSOLO_MM      = config["quant"].get("starsolo", {}).get("mm", "Unique")
TRIMMER          = config["quant"]["starsolo"]["trimmer"]
PREP             = config["quant"]["starsolo"]["preprocessor"]
if TRIMMER == 'starsolo':
    TRIMMER = ''

# splitpipe defaults (from libprep config)
CHEM = config["quant"]["chemistry"]
KIT =  config["quant"]["kit"].lower()


def _tso_window_args(tso: str, kmax: int = 15) -> str:
    # Emit cutadapt -g TSO{k}=^N{0..k}TSO for ED<=2-in-window style starts
    return " ".join([f"-g TSO{k}=^{'N'*k}{tso}" for k in range(kmax + 1)])


def star_extra_args(config):
    """
    Build STAR/STARsolo args (no read-order logic).

    Modes (config['quant']['starsolo']['trim']):
      - 'cutadapt' : trimming is upstream; STAR gets no adapter flags.
      - 'starsolo'      : STAR trims in CellRanger4 mode with ONLY 5' adapter (Parse TSO).
    """
    q  = config["quant"]
    ss = q.get("starsolo", {})

    kit = q["kit"].lower() # 'wt' | 'wt_mega' | 'wt_mini'
    preprocessor = ss.get("preprocessor", "").lower() # 'none' | 'rt_merge' | 'splitcode' 
    trimmer = ss.get("trimmer", "").lower() # 'none' | 'cutadapt' | 'starsolo'
    use_velo   = ss.get("use_velo", False)
    feature    = ss.get("feature_count", "GeneFull_Ex50pAS") # 'Gene', 'GeneFull_Ex50pAS'
    multimaps  = ss.get("multi_mappers", "Unique")  # 'Unique', 'EM'
    output_bam = ss.get("output_bam", False)

    # expected cells (kit-level, with optional override)
    n_by_kit = {"wt": 100_000, "wt_mega": 1_000_000, "wt_mini": 20_000}
    n_expected = int(q.get("n_expected_cells", n_by_kit[kit]))
    n_expected = int(n_expected/(len(SUBLIBS) * 1.1))

    args = ["--genomeLoad", "LoadAndKeep",
            "--soloCellReadStats", "Standard",
            "--soloBarcodeReadLength", "0",
            "--outSAMtype", "None",
            "--soloStrand", "Unstranded",
            "--soloCellFilter", "None"]

    # --soloFeatures: always include Gene; add feature if not 'Gene'; optionally add Velocyto
    feat_tokens = ["Gene"]
    if feature and feature != "Gene":
        # STAR expects tokens like 'GeneFull_Ex50pAS' as a single token
        feat_tokens.append(feature)
    args += ["--soloFeatures", *feat_tokens]
    if use_velo:
        args.append("Velocyto")
        
    if multimaps:
        args += ["--soloMultiMappers", str(multimaps)]

    if output_bam:
        args += ["--outSAMtype",  "BAM"]
        
    # format {preprocessor}
    mode = f"{preprocessor}_{trimmer}"
    if mode == "_":
        args += ["--soloCBmatchWLtype", "EditDist_2"]
    elif mode == "_cutadapt":
        args += ["--soloCBmatchWLtype", "EditDist_2"]
    elif mode == "_starsolo":
        args += ["--soloCBmatchWLtype", "EditDist_2"]
        args += ["--clipAdapterType", "CellRanger4",
                 "--clip5pAdapterSeq", PRE_TSO_SEQ]
    elif mode == "rt_merge_cutadapt":
        args += ["--soloCBmatchWLtype", "EditDist_2"]
    elif mode == "rt_merge_starsolo":
        args += ["--soloCBmatchWLtype", "EditDist_2"]
        args += ["--clipAdapterType", "CellRanger4",
                 "--clip5pAdapterSeq", PRE_TSO_SEQ]
    elif mode == "splitcode_starsolo":
        args += ["--soloCBmatchWLtype", "Exact"]
        args += ["--clipAdapterType", "CellRanger4",
                 "--clip5pAdapterSeq", PRE_TSO_SEQ]
    elif mode == "splitcode_cutadapt":
        args += ["--soloCBmatchWLtype", "1MM"]
    
    else:
        ValueError(f"Unsupported trim mode for {preprocessor} preprocessor: {trimmer}")

    return args


rule parsebio_ext:
    output:
        join(EXT_DIR, 'parsebio', '{name}_{kit}_{chem}.txt')
    log:
        join(EXT_DIR, 'parsebio', 'logs', '{name}_{kit}_{chem}.log')
    params:
        url = join(config.get('winecellar', {}).get('url', ''), 'parsebio', '{name}_{kit}_{chem}.txt'),
        proxy = config.get('proxy', {}).get('wget', '')
    shell:
        """
        wget {params.proxy} {params.url} -O- > {output}
        echo "Parse Biosciences {wildcards.name},NA,{params.url},`date -I`" > {log}
        """

        
rule parsebio_whitelists:
    input:
        barcodes = "/mnt/archive/ext_cache/ext/parsebio/barcodes/bc_data_v1.csv"
    params:
        script     = src_gcf("scripts/gen_whitelists.py"),
        barcodes_dir = "/mnt/archive/ext_cache/ext/parsebio/barcodes",
        kit        = KIT,
        chemistry  = CHEM,
        trimmer    = config.get("quant",{}).get("starsolo",{}).get("trim","starsolo"),
        outdir     = "data/tmp/singlecell/whitelists/",
    output:
        r1_R       = "data/tmp/singlecell/whitelists/r1_R.txt",
        r1_T       = "data/tmp/singlecell/whitelists/r1_T.txt",
        r1         = "data/tmp/singlecell/whitelists/r1.txt",
        r2         = "data/tmp/singlecell/whitelists/r2.txt",
        r3         = "data/tmp/singlecell/whitelists/r3.txt",
        r1_wm      = "data/tmp/singlecell/whitelists/r1_wellmap.txt",
        r2_wm      = "data/tmp/singlecell/whitelists/r2_wellmap.txt",
        r3_wm      = "data/tmp/singlecell/whitelists/r3_wellmap.txt",  
    shell:
        "python {params.script} "
        "--kit {params.kit} "
        "--chem {params.chemistry} "
        "--barcodes-dir {params.barcodes_dir} "
        "--outdir {params.outdir} "

# Reformat (existing) → config.txt
rule parsebio_splitcode_config_reformat:
    input:
        r1_R = "data/tmp/singlecell/whitelists/r1_R.txt",
        r1_T = "data/tmp/singlecell/whitelists/r1_T.txt",
        r2   = "data/tmp/singlecell/whitelists/r2.txt",
        r3   = "data/tmp/singlecell/whitelists/r3.txt",
    params:
        script    = src_gcf("scripts/splitcode_config.py"),
        chemistry = CHEM,
        outdir    = "data/tmp/singlecell/fastq/splitcode/",
        read_idx  = 1,     # R2 in paired runs
    output:
        config = "data/tmp/singlecell/fastq/splitcode/config.txt",
    shell:
        """
        python {params.script} \
          --mode reformat \
          --chem {params.chemistry} \
          --outdir {params.outdir} \
          --read-index {params.read_idx} \
          --r1-T {input.r1_T} \
          --r1-R {input.r1_R} \
          --r2 {input.r2} \
          --r3 {input.r3} \
        """

# Convert R→T (new) → config_rt.txt
rule parsebio_splitcode_config_rt:
    input:
        r1_R = "data/tmp/singlecell/whitelists/r1_R.txt",
        r1_T = "data/tmp/singlecell/whitelists/r1_T.txt",
    params:
        script    = src_gcf("scripts/splitcode_config.py"),
        chemistry = CHEM,
        outdir    = "data/tmp/singlecell/fastq/rt_merge",
        read_idx  = 1,   # keep consistent; override to 0 if running R2-only
        dist      = 1,
    output:
        config_rt = "data/tmp/singlecell/fastq/rt_merge/config.txt",
    shell:
        """
        python {params.script} \
          --mode rt-convert \
          --chem {params.chemistry} \
          --outdir {params.outdir} \
          --read-index {params.read_idx} \
          --rt-distance {params.dist} \
          --r1-R {input.r1_R} \
          --r1-T {input.r1_T}
        """


def tso_window_args(tso_seq, kmax=15):
    return " ".join([f"-g TSO{k}=^{'N'*k}{tso_seq}" for k in range(0, kmax+1)])


rule parsebio_fastq_raw:
    input:
        unpack(get_raw_fastq)
    output:
        R1 = "data/tmp/singlecell/fastq/{sublib}_R1.fastq.gz",
        R2 = "data/tmp/singlecell/fastq/{sublib}_R2.fastq.gz",
    wildcard_constraints:
        sublib = r"[^/]+"
    shell:
        "ln -srnf {input.R1} {output.R1};  ln -srnf {input.R2} {output.R2} "

rule parsebio_fastq_rt_merge:
    input:
        R1 = "data/tmp/singlecell/fastq/{sublib}_R1.fastq.gz",
        R2 = "data/tmp/singlecell/fastq/{sublib}_R2.fastq.gz",
        conf = "data/tmp/singlecell/fastq/rt_merge/config.txt"
    output:
        R1 = temp("data/tmp/singlecell/fastq/rt_merge/{sublib}_R1.fastq.gz"),
        R2 = temp("data/tmp/singlecell/fastq/rt_merge/{sublib}_R2.fastq.gz")
    wildcard_constraints:
        sublib = r"[^/]+"
    log:
        "data/tmp/singlecell/fastq/rt_merge/{sublib}.log" 
    threads:
        4
    group:
        groupname
    container:
        'docker://' + config['docker']['star']
    shell:
        "splitcode -c {input.conf} --summary {log} -t {threads} --nFastqs 2 -o {output.R1},{output.R2} {input.R1} {input.R2}"

rule parsebio_fastq_splitcode:
    input:
        R1 = "data/tmp/singlecell/fastq/{sublib}_R1.fastq.gz",
        R2 = "data/tmp/singlecell/fastq/{sublib}_R2.fastq.gz",
        conf = "data/tmp/singlecell/fastq/splitcode/config.txt"
    output:
        R1 = temp("data/tmp/singlecell/fastq/splitcode/{sublib}_R1.fastq.gz"),
        R2 = temp("data/tmp/singlecell/fastq/splitcode/{sublib}_R2.fastq.gz")
    wildcard_constraints:
        sublib = r"[^/]+"
    log:
        "data/tmp/singlecell/fastq/splitcode/{sublib}.log" 
    threads:
        4
    group:
        groupname
    container:
        'docker://' + config['docker']['star']
    shell:
        "splitcode -c {input.conf} --summary {log} -t {threads} --nFastqs 2 -o {output.R1},{output.R2} {input.R1} {input.R2}"
        

rule parsebio_fastq_trim_cutadapt:
    input:
        R1 = join("data/tmp/singlecell/fastq", PREP, "{sublib}_R1.fastq.gz"),
        R2 = join("data/tmp/singlecell/fastq", PREP, "{sublib}_R2.fastq.gz")
    output:
        R1 = temp(join("data/tmp/singlecell/fastq", PREP, "cutadapt", "{sublib}_R1.fastq.gz")),
        R2 = temp(join("data/tmp/singlecell/fastq", PREP, "cutadapt", "{sublib}_R2.fastq.gz")),
        json = "data/tmp/singlecell/fastq/cutadapt/{sublib}.log.json"
    wildcard_constraints:
        sublib = r"[^/]+"
    params:
        tso    = PRE_TSO_SEQ,
        tso_args = tso_window_args(PRE_TSO_SEQ),
        l12    = LINKER_RC,
        minlen_R1 = config["quant"]["starsolo"].get("min_r1_len", 25),
        minlen_R2 = config['read_geometry'][1],
        dirname = "data/tmp/singlecell/fastq/cutadapt/{sublib}"
    threads:
        12
    group:
        groupname
    container:
        'docker://' + config['docker']['star']
    shell:
        'cutadapt '
        '{params.tso_args} '
        '-a L12={params.l12} '
        '-a polyA=A{{15}}$ '
        '-n 2 --no-indels -e 0.15 -O 8 '
        '-m {params.minlen_R1}:{params.minlen_R2} '
        '--report=minimal --json {output.json} '
        '-j {threads} '
        '-o {output.R1} '
        '-p {output.R2} '
        ' {input.R1} '
        ' {input.R2} '


rule parsebio_starsolo_quant:
    input:
        R1 = join("data/tmp/singlecell/fastq", PREP, TRIMMER, "{sublib}_R1.fastq.gz"),
        R2 = join("data/tmp/singlecell/fastq", PREP, TRIMMER, "{sublib}_R2.fastq.gz"),
        genome = join(REF_DIR, 'index', 'genome', 'splitpipe', 'SA'),
        wl_1 = "data/tmp/singlecell/whitelists/r1.txt",
        wl_2 = "data/tmp/singlecell/whitelists/r2.txt",
        wl_3 = "data/tmp/singlecell/whitelists/r3.txt",
        #conf = "data/tmp/singlecell/fastq/splitcode/config_rt.txt"
    threads:
        48
    params:
        outdir = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}') + '/',
        genome_dir = join(REF_DIR, 'index', 'genome', 'splitpipe'),
        soloCBposition = config["quant"]["starsolo"]["soloCBposition"],
        soloUMIposition = config["quant"]["starsolo"]["soloUMIposition"],
        #pipe_block = preprocessor_call(config),
        #star_args = build_star_args(config)
    output:
        #mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'matrix.mtx'),
        #barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'barcodes.tsv'),
        #genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'features.tsv'),
        raw_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'matrix.mtx'),
        raw_barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'barcodes.tsv'),
        raw_genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'features.tsv'),
        #spliced_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'spliced.mtx'),
        #unspliced_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'unspliced.mtx'),
        #velo_barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'barcodes.tsv'),
        #velo_features = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'features.tsv'),
        #raw_mtx_em = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'UniqueAndMult-EM.mtx'),
        #bam = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Aligned.sortedByCoord.out.bam'),
        gene_stats = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'Features.stats'),
        barcode_stats = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'CellReads.stats'),
        summary_stats = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'Summary.csv'),
        #barcode_rank = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'UMIperCellSorted.txt'),
    container:
        'docker://' + config['docker']['star']
    benchmark:
        'benchmarks/parsebio_starsolo/{sublib}-starsolo.txt'
    log:
        star = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Log.final.out'),
        barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Barcodes.stats'),
        umi_cell = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'UMIperCellSorted.txt'),
        splitcode = "logs/{sublib}_splitcode.txt",
        cutadapt = "logs/{sublib}_cutadapt.txt"
    shell:
        'STAR --soloType CB_UMI_Complex '
        '--readFilesIn {input.R1} {input.R2} '
        '--readFilesCommand zcat '
        '--soloCBwhitelist {input.wl_3} {input.wl_2} {input.wl_1} '
        '--genomeDir {params.genome_dir} '
        '--outFileNamePrefix {params.outdir} '
        '--soloCBposition {params.soloCBposition} '
        '--soloUMIposition {params.soloUMIposition} '
        '--runThreadN {threads} '
        + 
        ' '.join(star_extra_args(config))


rule parsebio_starsolo_filtered:
    input:
        raw_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'matrix.mtx'),
        raw_barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'barcodes.tsv'),
        raw_genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'features.tsv')
    output:
        mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'matrix.mtx'),
        barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'barcodes.tsv'),
        genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'features.tsv')
    params:
        script = src_gcf("scripts/parsebio_barcode_rank.py"),
        n_expected_cells = 50_000
    container:
        'docker://' + config['docker']['default']
    shell:
        'python {params.script} '
        '--input-mtx {input.raw_mtx} '
        '--output-mtx {output.mtx} '


rule parsebio_starsolo_mtx_v2_fix:
    input:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, '{dge_type}', 'features.tsv')
    output:
        temp(join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, '{dge_type}', 'genes.tsv'))
    shell:
        "cp {input} {output}"


rule parsebio_starsolo_clean_shmem:
    input:
        expand(rules.parsebio_starsolo_quant.output, sublib=SUBLIBS)
    params:
        genome_dir = rules.parsebio_starsolo_quant.params.genome_dir
    output:
        temp(touch(join(QUANT_INTERIM, 'parsebio_starsolo', '.starsolo.mem.cleaned')))
    shadow:
        'minimal'
    container:
        'docker://' + config['docker']['star']
    shell:
        'STAR --genomeDir {params.genome_dir} --genomeLoad Remove || echo "no shared mem"'


rule parsebio_starsolo_barcode_sublib:
    input:
        r1_wm      = "data/tmp/singlecell/whitelists/r1_wellmap.txt",
        r2_wm      = "data/tmp/singlecell/whitelists/r2_wellmap.txt",
        r3_wm      = "data/tmp/singlecell/whitelists/r3_wellmap.txt",
        barcodes = join(QUANT_INTERIM, "parsebio_starsolo", "{sublib}", "Solo.out", STARSOLO_FEATURE, "raw", "barcodes.tsv")
    output:
        info = join(QUANT_INTERIM, "parsebio_starsolo", "{sublib}", "barcode_info.tsv")
    params:
        script = src_gcf("scripts/parsebio_starsolo_sampleinfo.py"),
        config = workflow.configfiles[0],
        sublib = "{sublib}"
    container:
        'docker://' + config['docker']['default']
    shell:
        "python {params.script} "
        "--r1-wellmap {input.r1_wm} "
        "--r2-wellmap {input.r2_wm} "
        "--r3-wellmap {input.r3_wm} "
        "--configfile {params.config} "
        "--sublib {params.sublib} "
        "--barcodes {input.barcodes} "
        "--order r3_r2_r1 "
        "--output {output.info} "

rule parsebio_starsolo_barcode_info:
    input:
        expand(join(QUANT_INTERIM, "parsebio_starsolo", "{sublib}", "barcode_info.tsv"), sublib=SUBLIBS)
    output:
        info = join(QUANT_INTERIM, "parsebio_starsolo", "barcode_info.tsv")
    container:
        'docker://' + config['docker']['default']
    shell:
        r"""awk 'FNR==1 && NR!=1 {{next}}; 1' {input} > {output.info}"""

        
rule parsebio_starsolo_scanpy_pp_ipynb:
    input:
        join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', 'scanpy', '{aggr_id}_filtered.h5ad')
    output:
        preprocessed = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', 'scanpy', '{aggr_id}_preprocessed.h5ad'),
    log:
        notebook = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    threads:
        24
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    notebook:
        'scripts/parsebio_starsolo_preprocess.py.ipynb'


rule parsebio_starsolo_scanpy_pp_ipynb_html:
    input:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_preprocessed.h5ad')
    output:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', 'notebooks', '{aggr_id}_pp.html')
    params:
        notebook = join(QUANT_INTERIM, 'aggregate', '{method}' 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')

    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.notebook} ' 
