#-*- mode:snakemake -*-
"""

"""

import re
import tempfile

include:
    join(GCFDB_DIR, 'parsebio.db')

# Parse adapters
PRE_TSO_SEQ      = config["quant"].get("tso", "AACGCAGAGTGAATGGG")
LINKER_RC        = config["quant"].get("l21_rc", "AACGCAGAGTGAATGGG")
STARSOLO_TRIMMER = STARSOLO_CONFIG.get("trimmer", "skip")
TRIMMER          = "" if STARSOLO_TRIMMER == "starsolo" else STARSOLO_TRIMMER
PREP             = STARSOLO_CONFIG.get("preprocessor", "skip")
if STARSOLO_OUTPUT_BAM:
    STARSOLO_PARSEBIO_BAM_OUTPUT = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Aligned.sortedByCoord.out.bam')
else:
    STARSOLO_PARSEBIO_BAM_OUTPUT = []

# splitpipe defaults (from libprep config)
CHEM = config["quant"]["chemistry"]
KIT =  config["quant"]["kit"].lower()
n_by_kit = {"wt": 100_000, "wt_mega": 1_000_000, "wt_mini": 20_000}
n_expected = int(config["quant"].get("n_expected_cells", n_by_kit[KIT]))
try:
    n_sublibs = len(SUBLIBS)  # noqa: F821
except NameError:
    n_sublibs = 1
N_EXPECTED_CELLS = max(1, int(n_expected / (n_sublibs * 1.1)))
max_cells_by_kit = {"wt_mini": 30_000,
                    "wt": 30_000,
                    "wt_mega": 150_000,
                    "wt_mega_384": 300_000,
                    "wt_penta": 300_000,
                    "wt_penta_384": 500_000,
                    "custom": 150_000,
}
BARCODE_RANK_MAX_CELLS = int(config["quant"].get("barcode_rank_max_cells", max_cells_by_kit[KIT]))
STARSOLO_PARSEBIO_MITO_NAMES = STARSOLO_MITO_NAMES.copy()

db_conf = config["db"][config["db"]["reference_db"]]
assembly = db_conf.get("assembly", config["db"].get("assembly"))
if not assembly:
    raise ValueError("No 'assembly' found for selected reference_db in config")

def _parsebio_genome_name(name: str) -> str:
    """
    Reproduce split-pipe utils.clean_vname_chars(name, to_dash=True).
    Allowed: A–Z, a–z, 0–9, and '-'.
    All other chars (including '.' and '_') -> '-'.
    Consecutive '-' collapsed, leading/trailing '-' stripped.
    """
    if not name:
        raise ValueError("Empty genome assembly name")

    # Replace all non [A-Za-z0-9-] with '-'
    cleaned = re.sub(r"[^A-Za-z0-9-]+", "-", name)
    # Collapse multiple '-' into one
    cleaned = re.sub(r"-+", "-", cleaned)
    # Strip leading/trailing '-'
    cleaned = cleaned.strip("-")

    if not cleaned:
        raise ValueError(f"Assembly name '{name}' cleans to empty under ParseBio rules")

    return cleaned

pb_assembly_name = _parsebio_genome_name(assembly)
STARSOLO_PARSEBIO_MITO_NAMES += [f"{pb_assembly_name}_{name}" for name in STARSOLO_MITO_NAMES]

STARSOLO_PARSEBIO_ARGS = STARSOLO_COMMON_ARGS + [
    "--soloType", "CB_UMI_Complex",
    "--readFilesCommand", "zcat",
    "--soloBarcodeReadLength", "0",
    "--soloStrand", "Unstranded",
    "--soloCellFilter", "None",
    "--outSAMmultNmax", "3",
    "--soloUMIdedup", STARSOLO_PARSEBIO_UMI_DEDUP,
    "--soloUMIfiltering", STARSOLO_PARSEBIO_UMI_FILTERING,
    "--genomeChrSetMitochondrial", *STARSOLO_PARSEBIO_MITO_NAMES,
]

mode = f"{PREP}_{STARSOLO_TRIMMER}"
if mode in {
    "skip_skip",
    "skip_cutadapt",
    "skip_starsolo",
    "rt_merge_cutadapt",
    "rt_merge_starsolo",
    "error_correct_bc1_cutadapt",
    "error_correct_bc1_starsolo",
}:
    STARSOLO_PARSEBIO_ARGS += ["--soloCBmatchWLtype", "EditDist_2"]
    if mode.endswith("_starsolo"):
        STARSOLO_PARSEBIO_ARGS += ["--clipAdapterType", "CellRanger4", "--clip5pAdapterSeq", PRE_TSO_SEQ]
elif mode == "splitcode_starsolo":
    STARSOLO_PARSEBIO_ARGS += [
        "--soloCBmatchWLtype", "Exact",
        "--clipAdapterType", "CellRanger4",
        "--clip5pAdapterSeq", PRE_TSO_SEQ,
    ]
elif mode == "splitcode_cutadapt":
    STARSOLO_PARSEBIO_ARGS += ["--soloCBmatchWLtype", "1MM"]
else:
    raise ValueError(f"Unsupported trim mode: preprocessor={PREP}, trimmer={STARSOLO_TRIMMER}")

STARSOLO_PARSEBIO_ARGS = " ".join(STARSOLO_PARSEBIO_ARGS)

if VELO_OUTPUT:
    STARSOLO_PARSEBIO_VELOCYTO_RAW_BARCODES = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'barcodes.tsv')
    STARSOLO_PARSEBIO_VELOCYTO_RAW_FEATURES = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'features.tsv')
    STARSOLO_PARSEBIO_VELOCYTO_RAW_SPLICED = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'spliced.mtx')
    STARSOLO_PARSEBIO_VELOCYTO_RAW_UNSPLICED = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'unspliced.mtx')
    STARSOLO_PARSEBIO_VELOCYTO_RAW_AMBIGUOUS = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Velocyto', 'raw', 'ambiguous.mtx')
else:
    STARSOLO_PARSEBIO_VELOCYTO_RAW_BARCODES = []
    STARSOLO_PARSEBIO_VELOCYTO_RAW_FEATURES = []
    STARSOLO_PARSEBIO_VELOCYTO_RAW_SPLICED = []
    STARSOLO_PARSEBIO_VELOCYTO_RAW_UNSPLICED = []
    STARSOLO_PARSEBIO_VELOCYTO_RAW_AMBIGUOUS = []



def _tso_window_args(tso: str, kmax: int = 15) -> str:
    # Emit cutadapt -g TSO{k}=^N{0..k}TSO for ED<=2-in-window style starts
    return " ".join([f"-g TSO{k}=^{'N'*k}{tso}" for k in range(kmax + 1)])


def star_build_preprocessor_string(config, pipe=False):
    q  = config["quant"]
    ss = q.get("starsolo", {})
    preprocessor = ss.get("preprocessor", "").lower()
    trimmer      = ss.get("trimmer", "").lower()

    tmpdir = tempfile.mkdtemp(prefix="snakemake-fifo-")
    atexit.register(lambda: shutil.rmtree(tmpdir, ignore_errors=True))
    tmp_pp_r1 = os.path.join(tmpdir, "pp_R1.fastq")
    tmp_pp_r2 = os.path.join(tmpdir, "pp_R2.fastq")
    if pipe:
        os.mkfifo(tmp_pp_r1)
        os.mkfifo(tmp_pp_r2)

    if preprocessor in  ["rt_merge", "splitcode", "error_correct_bc1"]:
        config_file = join(FILTER_INTERIM, "fastq", preprocessor, "config.txt")
        log_file = join(FILTER_INTERIM, "fastq", preprocessor, "log.txt")
        args = f"splitcode -c {config_file} --summary {log_file} -t {{threads}} --nFastqs 2 -o {tmp_pp_r1},{tmp_pp_r2} {{input.R1}} {{input.R2}};\n"
    elif preprocessor == "skip":
        args = f"zcat {{input.R1}} > {tmp_pp_r1}; \nzcat {{input.R2}} > {tmp_pp_r2}; \n"
    else:
        raise ValueError
    
    if trimmer == "cutadapt":
        tmp_trim_r1 = os.path.join(tmpdir, "trim_R1.fastq")
        tmp_trim_r2 = os.path.join(tmpdir, "trim_R2.fastq")
        if pipe:
            os.mkfifo(tmp_trim_r1)
            os.mkfifo(tmp_trim_r2)
        tso_args = _tso_window_args(PRE_TSO_SEQ)
        
        args += f'cutadapt {tso_args} -a L12={LINKER_RC} -a polyA=A{{15}}$ -n 2 --no-indels -e 0.15 -O 8 -m {minlen_R1}:{minlen_R2} --report=minimal --json {output.json} -j {{threads}} '
        '-o {output.R1} '
        '-p {output.R2} '
        ' {input.R1} '
        ' {input.R2} '

    elif trimmer == "skip":
        pass
    else:
        raise ValueError


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
        barcodes = parsebio_barcode_inputs(KIT, CHEM)
    params:
        script     = src_gcf("scripts/gen_whitelists.py"),
        barcodes_dir = PARSEBIO_BARCODE_DIR,
        kit        = KIT,
        chemistry  = CHEM,
        trimmer    = config.get("quant",{}).get("starsolo",{}).get("trim","starsolo"),
        outdir     = join(INTERIM_DIR, "singlecell", "whitelists"),
    output:
        r1_R       = join(INTERIM_DIR, "singlecell", "whitelists", "r1_R.txt"),
        r1_T       = join(INTERIM_DIR, "singlecell", "whitelists", "r1_T.txt"),
        r1         = join(INTERIM_DIR, "singlecell", "whitelists", "r1.txt"),
        r2         = join(INTERIM_DIR, "singlecell", "whitelists", "r2.txt"),
        r3         = join(INTERIM_DIR, "singlecell", "whitelists", "r3.txt"),
        r1_wm      = join(INTERIM_DIR, "singlecell", "whitelists", "r1_wellmap.txt"),
        r2_wm      = join(INTERIM_DIR, "singlecell", "whitelists", "r2_wellmap.txt"),
        r3_wm      = join(INTERIM_DIR, "singlecell", "whitelists", "r3_wellmap.txt"),  
    shell:
        "python {params.script} "
        "--kit {params.kit} "
        "--chem {params.chemistry} "
        "--barcodes-dir {params.barcodes_dir} "
        "--outdir {params.outdir} "

# Reformat (existing) → config.txt
rule parsebio_splitcode_config_reformat:
    input:
        r1_R = join(INTERIM_DIR, "singlecell", "whitelists", "r1_R.txt"),
        r1_T = join(INTERIM_DIR, "singlecell", "whitelists", "r1_T.txt"),
        r2   = join(INTERIM_DIR, "singlecell", "whitelists", "r2.txt"),
        r3   = join(INTERIM_DIR, "singlecell", "whitelists", "r3.txt"),
    params:
        script    = src_gcf("scripts/splitcode_config.py"),
        chemistry = CHEM,
        outdir    = join(FILTER_INTERIM, "fastq", "splitcode"),
        read_idx  = 1,     # R2 in paired runs
    output:
        config = join(FILTER_INTERIM, "fastq", "splitcode/config.txt"),
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

# Convert R→T
rule parsebio_splitcode_config_rt:
    input:
        r1_R = join(INTERIM_DIR, "singlecell", "whitelists", "r1_R.txt"),
        r1_T = join(INTERIM_DIR, "singlecell", "whitelists", "r1_T.txt"),
    params:
        script    = src_gcf("scripts/splitcode_config.py"),
        chemistry = CHEM,
        outdir    = join(FILTER_INTERIM, "fastq", "rt_merge"),
        read_idx  = 1,   # keep consistent; override to 0 if running R2-only
        dist      = 1,
    output:
        config_rt = join(FILTER_INTERIM, "fastq", "rt_merge", "config.txt")
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

rule parsebio_splitcode_error_correct_bc1:
    input:
        r1_R = join(INTERIM_DIR, "singlecell", "whitelists", "r1_R.txt"),
        r1_T = join(INTERIM_DIR, "singlecell", "whitelists", "r1_T.txt"),
    params:
        script    = src_gcf("scripts/splitcode_config.py"),
        chemistry = CHEM,
        outdir    = join(FILTER_INTERIM, "fastq", "error_correct_bc1"),
        read_idx  = 1,   # keep consistent; override to 0 if running R2-only
        dist      = 1,
    output:
        config_rt = join(FILTER_INTERIM, "fastq", "error_correct_bc1", "config.txt")
    shell:
        """
        python {params.script} \
          --mode error-correct-bc1 \
          --chem {params.chemistry} \
          --outdir {params.outdir} \
          --read-index {params.read_idx} \
          --rt-distance {params.dist} \
          --r1-R {input.r1_R} \
          --r1-T {input.r1_T}
        """
    

def tso_window_args(tso_seq, kmax=15):
    return " ".join([f"-g TSO{k}=^{'N'*k}{tso_seq}" for k in range(0, kmax+1)])



rule parsebio_fastq_rt_merge:
    input:
        R1 = join(FILTER_INTERIM, "fastq", "{sublib}_R1.fastq.gz"),
        R2 = join(FILTER_INTERIM, "fastq", "{sublib}_R2.fastq.gz"),
        conf = join(FILTER_INTERIM, "fastq", "rt_merge", "config.txt")
    output:
        R1 = temp(join(FILTER_INTERIM, "fastq", "rt_merge", "{sublib}_R1.fastq.gz")),
        R2 = temp(join(FILTER_INTERIM, "fastq", "rt_merge", "{sublib}_R2.fastq.gz"))
    wildcard_constraints:
        sublib = r"[^/]+"
    log:
        join(FILTER_INTERIM, "fastq", "rt_merge", "{sublib}.log") 
    threads:
        12
    container:
        'docker://' + config['docker']['star']
    shell:
        "splitcode -c {input.conf} --summary {log} -t {threads} --nFastqs 2 -o {output.R1},{output.R2} {input.R1} {input.R2}"

        
rule parsebio_fastq_error_correct_bc1:
    input:
        R1 = join(FILTER_INTERIM, "fastq", "{sublib}_R1.fastq.gz"),
        R2 = join(FILTER_INTERIM, "fastq", "{sublib}_R2.fastq.gz"),
        conf = join(FILTER_INTERIM, "fastq", "error_correct_bc1", "config.txt")
    output:
        R1 = temp(join(FILTER_INTERIM, "fastq", "error_correct_bc1", "{sublib}_R1.fastq.gz")),
        R2 = temp(join(FILTER_INTERIM, "fastq", "error_correct_bc1", "{sublib}_R2.fastq.gz"))
    wildcard_constraints:
        sublib = r"[^/]+"
    log:
        join(FILTER_INTERIM, "fastq", "error_correct_bc1", "{sublib}.log") 
    threads:
        12
    container:
        'docker://' + config['docker']['star']
    shell:
        "splitcode -c {input.conf} --summary {log} -t {threads} --nFastqs 2 -o {output.R1},{output.R2} {input.R1} {input.R2}"

rule parsebio_fastq_splitcode:
    input:
        R1 = join(FILTER_INTERIM, "fastq", "{sublib}_R1.fastq.gz"),
        R2 = join(FILTER_INTERIM, "fastq", "{sublib}_R2.fastq.gz"),
        conf = join(FILTER_INTERIM, "fastq", "splitcode", "config.txt")
    output:
        R1 = temp(join(FILTER_INTERIM, "fastq", "splitcode", "{sublib}_R1.fastq.gz")),
        R2 = temp(join(FILTER_INTERIM, "fastq", "splitcode", "{sublib}_R2.fastq.gz"))
    wildcard_constraints:
        sublib = r"[^/]+"
    log:
        join(FILTER_INTERIM, "fastq", "splitcode", "{sublib}.log") 
    threads:
        4
    container:
        'docker://' + config['docker']['star']
    shell:
        "splitcode -c {input.conf} --summary {log} -t {threads} --nFastqs 2 -o {output.R1},{output.R2} {input.R1} {input.R2}"
        

rule parsebio_fastq_preprocessor_skip:
    input:
        R1 = join(FILTER_INTERIM, "fastq", "{sublib}_R1.fastq.gz"),
        R2 = join(FILTER_INTERIM, "fastq", "{sublib}_R2.fastq.gz"),
    output:
        R1 = temp(join(FILTER_INTERIM, "fastq", "skip", "{sublib}_R1.fastq.gz")),
        R2 = temp(join(FILTER_INTERIM, "fastq", "skip", "{sublib}_R2.fastq.gz")),
    wildcard_constraints:
        sublib = r"[^/]+"
    shell:
        "ln -srnf {input.R1} {output.R1};  ln -srnf {input.R2} {output.R2} "


rule parsebio_fastq_trimmer_skip:
    input:
        R1 = join(FILTER_INTERIM, "fastq", PREP, "{sublib}_R1.fastq.gz"),
        R2 = join(FILTER_INTERIM, "fastq", PREP, "{sublib}_R2.fastq.gz")
    output:
        R1 = temp(join(FILTER_INTERIM, "fastq", PREP, "skip", "{sublib}_R1.fastq.gz")),
        R2 = temp(join(FILTER_INTERIM, "fastq", PREP, "skip", "{sublib}_R2.fastq.gz")),
    wildcard_constraints:
        sublib = r"[^/]+"
    shell:
        "ln -srnf {input.R1} {output.R1};  ln -srnf {input.R2} {output.R2} "


rule parsebio_fastq_trim_cutadapt:
    input:
        R1 = join(FILTER_INTERIM, "fastq", PREP, "{sublib}_R1.fastq.gz"),
        R2 = join(FILTER_INTERIM, "fastq", PREP, "{sublib}_R2.fastq.gz")
    output:
        R1 = temp(join(FILTER_INTERIM, "fastq", PREP, "cutadapt", "{sublib}_R1.fastq.gz")),
        R2 = temp(join(FILTER_INTERIM, "fastq", PREP, "cutadapt", "{sublib}_R2.fastq.gz")),
        json = join(FILTER_INTERIM, "fastq", PREP, "cutadapt", "{sublib}.log.json")
    wildcard_constraints:
        sublib = r"[^/]+"
    params:
        tso    = PRE_TSO_SEQ,
        tso_args = tso_window_args(PRE_TSO_SEQ),
        l12    = LINKER_RC,
        minlen_R1 = config["quant"]["starsolo"].get("min_r1_len", 25),
        minlen_R2 = config['read_geometry'][1],
        dirname = join(FILTER_INTERIM, "fastq", "cutadapt", "{sublib}")
    threads:
        12
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


def get_parsebio_starsolo_genome():
    return join(REF_DIR, 'index', 'genome', 'splitpipe', 'SA')

def get_parsebio_starsolo_config():
    q  = config["quant"]
    ss = q.get("starsolo", {})
    preprocessor = ss.get("preprocessor", "").lower()
    return join(FILTER_INTERIM, "fastq", f"{preprocessor}", "config.txt")


rule parsebio_starsolo_quant:
    input:
        R1 = join(FILTER_INTERIM, "fastq", PREP, TRIMMER, "{sublib}_R1.fastq.gz"),
        R2 = join(FILTER_INTERIM, "fastq", PREP, TRIMMER, "{sublib}_R2.fastq.gz"),
        wl_1 = join(INTERIM_DIR, "singlecell", "whitelists", "r1.txt"),
        wl_2 = join(INTERIM_DIR, "singlecell", "whitelists", "r2.txt"),
        wl_3 = join(INTERIM_DIR, "singlecell", "whitelists", "r3.txt"),
        genome = get_parsebio_starsolo_genome()
    threads:
        48
    wildcard_constraints:
        sublib = r"[^/]+"
    params:
        outdir = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}') + '/',
        genome_dir = lambda wildcards, input: os.path.dirname(input.genome),
        soloCBposition = STARSOLO_CONFIG["soloCBposition"],
        soloUMIposition = STARSOLO_CONFIG["soloUMIposition"],
        starsolo_args = STARSOLO_PARSEBIO_ARGS
    output:
        raw_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', STARSOLO_MTX),
        raw_barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'barcodes.tsv'),
        raw_genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'features.tsv'),
        gene_stats = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'Features.stats'),
        gene_summary = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'Summary.csv'),
        cell_reads = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'CellReads.stats'),
        bam = STARSOLO_PARSEBIO_BAM_OUTPUT,
        spliced_mtx = STARSOLO_PARSEBIO_VELOCYTO_RAW_SPLICED,
        unspliced_mtx = STARSOLO_PARSEBIO_VELOCYTO_RAW_UNSPLICED,
        ambiguous_mtx = STARSOLO_PARSEBIO_VELOCYTO_RAW_AMBIGUOUS,
        velo_barcodes = STARSOLO_PARSEBIO_VELOCYTO_RAW_BARCODES,
        velo_features = STARSOLO_PARSEBIO_VELOCYTO_RAW_FEATURES,
    container:
        'docker://' + config['docker']['star']
    benchmark:
        'benchmarks/parsebio_starsolo/{sublib}-starsolo.txt'
    log:
        star = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Log.final.out'),
        barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', 'Barcodes.stats'),
        umi_cell = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'UMIperCellSorted.txt')
    shell:
        'STAR '
        '--readFilesIn {input.R1} {input.R2} '
        '--soloCBwhitelist {input.wl_3} {input.wl_2} {input.wl_1} '
        '--genomeDir {params.genome_dir} '
        '--outFileNamePrefix {params.outdir} '
        '--soloCBposition {params.soloCBposition} '
        '--soloUMIposition {params.soloUMIposition} '
        '--runThreadN {threads} '
        '{params.starsolo_args} '

rule parsebio_starsolo_filtered:
    input:
        raw_mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', STARSOLO_MTX),
        raw_barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'barcodes.tsv'),
        raw_genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'raw', 'features.tsv'),
        bc_info = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'barcode_info.tsv')
    output:
        mtx = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'matrix.mtx'),
        barcodes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'barcodes.tsv'),
        genes = join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, 'filtered', 'features.tsv')
    params:
        script = src_gcf('scripts/parsebio_barcode_rank.py'),
        n_expected_cells = N_EXPECTED_CELLS,
        max_cells = BARCODE_RANK_MAX_CELLS,
        cnt_scale_fac = float(config['quant'].get('barcode_rank_cnt_scale_fac', 0.70)),
        method = config['quant'].get('barcode_rank_method', 'C'), 
    container:
        'docker://' + config['docker']['default']
    threads:
        8
    shell:
        'python {params.script} '
        '--n-expected-cells {params.n_expected_cells} '
        '--max-cells {params.max_cells} '
        '--cnt-scale-fac {params.cnt_scale_fac} '
        '--method {params.method} '
        '--input-mtx {input.raw_mtx} '
        '--output-mtx {output.mtx} '


def parsebio_starsolo_rt_inputs(wc):
    sublibs = AGGR_IDS[wc.aggr_id]
    if CB_OUTPUT:
        inputs = [join(QUANT_INTERIM, wc.method, s, 'cellbender', f'{s}_filtered.h5') for s in sublibs]
    else:
        inputs = [get_filtered_mtx(SimpleNamespace(method=wc.method, sublib=s, sample=s))['mtx'] for s in sublibs]
    return {
        'inputs': inputs,
        'feature_info': get_feature_info_list(wc),
        'barcode_info': get_barcode_info_list(wc),
    }


rule parsebio_starsolo_scanpy_rt_filtered:
    wildcard_constraints:
        method = 'parsebio_starsolo'
    input:
        unpack(parsebio_starsolo_rt_inputs)
    params:
        script    = src_gcf('scripts/convert_scanpy.py'),
        bc_type   = lambda wc: BC_RENAME[wc.method],
        enable_cb = '--enable-cellbender' if CB_OUTPUT  else ''
    output:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'cellbender', 'scanpy', '{aggr_id}_rt.h5ad') if CB_OUTPUT else join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', '{aggr_id}_rt.h5ad')
    container:
        'docker://' + config['docker']['scanpy'],
    threads: 48
    log:
        join(QUANT_INTERIM, 'aggregate', '{method}', 'scanpy', 'logs', '{aggr_id}.log'),
    shell:
        'python {params.script} ' 
        '{input.inputs} '
        '--feature-info {input.feature_info} '
        '--barcode-info {input.barcode_info} '
        '--barcode-rename {params.bc_type} '
        '-o {output} '
        '-f {wildcards.method} '
        '--log {log} '
        '{params.enable_cb} '
        '--verbose '


rule parsebio_starsolo_mtx_v2_fix:
    input:
        join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, '{dge_type}', 'features.tsv')
    output:
        temp(join(QUANT_INTERIM, 'parsebio_starsolo', '{sublib}', 'Solo.out', STARSOLO_FEATURE, '{dge_type}', 'genes.tsv'))
    shell:
        'cp {input} {output}'


rule parsebio_starsolo_clean_shmem:
    input:
        expand(rules.parsebio_starsolo_quant.output.raw_mtx, sublib=SUBLIBS)
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
        r1_wm      = join(INTERIM_DIR, "singlecell", "whitelists", "r1_wellmap.txt"),
        r2_wm      = join(INTERIM_DIR, "singlecell", "whitelists", "r2_wellmap.txt"),
        r3_wm      = join(INTERIM_DIR, "singlecell", "whitelists", "r3_wellmap.txt"),
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
        join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', 'scanpy', '{aggr_id}_preprocessed.h5ad')
    output:
        join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', 'scanpy', 'notebooks', '{aggr_id}_pp.html')
    params:
        notebook = join(QUANT_INTERIM, 'aggregate', 'parsebio_starsolo', 'scanpy', 'notebooks', '{aggr_id}_pp.ipynb')
    container:
        'docker://' + config['docker']['jupyter-scanpy']
    threads:
        1
    shell:
        'jupyter nbconvert --to html {params.notebook} '
