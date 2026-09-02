# allen_institute.smk — ABC downloads (wget + provenance) and derived palettes/maps
import re
from os.path import join
from urllib.parse import quote as _q

# -----------------------------
# Config shorthands
# -----------------------------
ABC = config['db']['allen_abc']
S3  = ABC.get('s3', {})
BUCKET = S3.get('bucket', 'allen-brain-cell-atlas').strip('/')
REGION = S3.get('region', 'us-west-2')
S3_HOST = f"https://{BUCKET}.s3.{REGION}.amazonaws.com"

WGET_PROXY = config.get('proxy', {}).get('wget', '')
DOCKER_DEF = 'docker://' + config['docker']['default']

# Output roots
EXT_ABC     = join(EXT_DIR, "allen-brain-cell-atlas")
ABC_MAPMY   = join(EXT_ABC, "mapmycells")
ABC_META    = join(EXT_ABC, "metadata")
ABC_DERIVED = join(EXT_ABC, "derived")

# -----------------------------
# URL helpers
# -----------------------------
def url_mapmycells(species: str, kind: str) -> str:
    """mapmycells/<product>/<release>/<file> -> full HTTPS URL"""
    s = ABC['mapmycells'][species]
    product = s['product']
    release = s['release']
    fname   = s['files'][kind]  # 'precomputed' or 'markers'
    key = f"mapmycells/{product}/{release}/{fname}"
    return f"{S3_HOST}/{_q(key, safe='/')}"

def url_taxonomy(file: str) -> str:
    """metadata/WMB-taxonomy/<release>/<file> -> full HTTPS URL"""
    rel = ABC['taxonomy']['release']
    key = f"metadata/WMB-taxonomy/{rel}/{file}"
    return f"{S3_HOST}/{_q(key, safe='/')}"

def url_metadata(product: str, release: str, file: str) -> str:
    """metadata/<PRODUCT>/<RELEASE>/<file> -> full HTTPS URL"""
    key = f"metadata/{product}/{release}/{file}"
    return f"{S3_HOST}/{_q(key, safe='/')}"

# -----------------------------
# MapMyCells: mouse (precomputed + markers)
# -----------------------------
rule mapmycells_precomputed_mouse:
    params:
        url     = url_mapmycells("mus_musculus", "precomputed"),
        release = ABC['mapmycells']['mus_musculus']['release'],
        name    = "ABC precomputed (mouse)",
        proxy   = WGET_PROXY
    output:
        pre_stats = join(ABC_MAPMY, "mus_musculus", "precomputed_stats.h5")
    log:
        join(ABC_MAPMY, "mus_musculus", "logs", "precomputed.prov.csv")
    shell:
        r'''
        wget {params.proxy} -O "{output.pre_stats}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''

rule mapmycells_markers_mouse:
    params:
        url     = url_mapmycells("mus_musculus", "markers"),
        release = ABC['mapmycells']['mus_musculus']['release'],
        name    = "ABC markers (mouse)",
        proxy   = WGET_PROXY
    output:
        markers = join(ABC_MAPMY, "mus_musculus", "markers.json")
    log:
        join(ABC_MAPMY, "mus_musculus", "logs", "markers.prov.csv")
    shell:
        r'''
        wget {params.proxy} -O "{output.markers}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''

# -----------------------------
# MapMyCells: human (precomputed + markers)
# -----------------------------
rule mapmycells_precomputed_human:
    params:
        url     = url_mapmycells("homo_sapiens", "precomputed"),
        release = ABC['mapmycells']['homo_sapiens']['release'],
        name    = "ABC precomputed (human)",
        proxy   = WGET_PROXY
    output:
        pre_stats = join(ABC_MAPMY, "homo_sapiens", "precomputed_stats.h5")
    log:
        join(ABC_MAPMY, "homo_sapiens", "logs", "precomputed.prov.csv")
    shell:
        r'''
        wget {params.proxy} -O "{output.pre_stats}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''

rule mapmycells_markers_human:
    params:
        url     = url_mapmycells("homo_sapiens", "markers"),
        release = ABC['mapmycells']['homo_sapiens']['release'],
        name    = "ABC markers (human)",
        proxy   = WGET_PROXY
    output:
        markers = join(ABC_MAPMY, "homo_sapiens", "markers.json")
    log:
        join(ABC_MAPMY, "homo_sapiens", "logs", "markers.prov.csv")
    shell:
        r'''
        wget {params.proxy} -O "{output.markers}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''

# -----------------------------
# Taxonomy (colors)
# -----------------------------
COLOR_RELEASE = ABC['taxonomy']['release']

rule abc_taxonomy_term_http:
    params:
        url     = url_taxonomy("cluster_annotation_term.csv"),
        release = COLOR_RELEASE,
        name    = "ABC taxonomy term table",
        proxy   = WGET_PROXY
    output:
        term_csv = join(ABC_META, "WMB-taxonomy", COLOR_RELEASE, "cluster_annotation_term.csv")
    log:
        join(ABC_META, "WMB-taxonomy", COLOR_RELEASE, "logs", "term.prov.csv")
    shell:
        r'''
        wget {params.proxy} -O "{output.term_csv}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''

# Build all-level colors (JSON + long TSV) — script runs in container
ABC_COLORS_LONG = join(ABC_DERIVED, "abc_colors.long.tsv")
ABC_COLORS_JSON = join(ABC_DERIVED, "abc_colors.json")

rule abc_build_colors:
    input:
        term_csv = rules.abc_taxonomy_term_http.output.term_csv
    output:
        long = ABC_COLORS_LONG,
        json = ABC_COLORS_JSON
    container:
        DOCKER_DEF
    shell:
        r'''
        python scripts/build_abc_colors_from_taxonomy.py \
          --term-csv "{input.term_csv}" \
          --out-long "{output.long}" \
          --out-json "{output.json}" \
          --levels neurotransmitter class subclass supertype cluster
        '''

# -----------------------------
# ROI metadata (dataset product) + roi_map.tsv
# -----------------------------
ROI = ABC.get('roi', {})
ROI_PRODUCT = ROI.get('product', ABC['mapmycells']['mus_musculus']['product'])
ROI_RELEASE = ROI.get('release', ABC['mapmycells']['mus_musculus']['release'])

rule abc_roi_metadata:
    params:
        url     = url_metadata(ROI_PRODUCT, ROI_RELEASE, "region_of_interest_metadata.csv"),
        release = ROI_RELEASE,
        name    = f"ABC ROI metadata ({ROI_PRODUCT})",
        proxy   = WGET_PROXY
    output:
        roi_csv = join(ABC_META, ROI_PRODUCT, ROI_RELEASE, "region_of_interest_metadata.csv")
    log:
        join(ABC_META, ROI_PRODUCT, ROI_RELEASE, "logs", "roi_metadata.prov.csv")
    shell:
        r'''
        wget {params.proxy} -O "{output.roi_csv}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''

rule abc_build_roi_map:
    input:
        roi = rules.abc_roi_metadata.output.roi_csv
    params:
        script = src_gcf("scripts/build_abc_roi_map.py")
    output:
        out = join(ABC_DERIVED, "roi_map.tsv")
    container:
        DOCKER_DEF
    shell:
        'python {params.script} '
         ' --roi {input.roi} '
         ' --out {output.out} '


# --- Neighbors ---------------------------------------------------------------
NEIGH = ABC.get('neighbors', {}) or {}
NEIGH_PRODUCT = NEIGH.get('product', ABC['mapmycells']['mus_musculus']['product'])
NEIGH_RELEASE = NEIGH.get('release', ROI_RELEASE)
NEIGH_FILES   = NEIGH.get('files') or ({'neighbors': NEIGH['file']} if NEIGH.get('file') else {})
NEIGH_FILE_PAT = "(" + "|".join(re.escape(v) for v in NEIGH_FILES.values()) + ")" if NEIGH_FILES else r"^$"
NEIGH_BASE_URL = f"{S3_HOST}/metadata/{NEIGH_PRODUCT}/{NEIGH_RELEASE}"

if NEIGH_FILES:
    rule abc_neighbors_file:
        wildcard_constraints:
            file = NEIGH_FILE_PAT
        params:
            proxy   = WGET_PROXY,
            release = NEIGH_RELEASE,
            base    = NEIGH_BASE_URL
        output:
            out = join(ABC_META, NEIGH_PRODUCT, NEIGH_RELEASE, "{file}")
        log:
            join(ABC_META, NEIGH_PRODUCT, NEIGH_RELEASE, "logs", "{file}.prov.csv")
        shell:
            r'''
            wget {params.proxy} -O "{output.out}" "{params.base}/{wildcards.file}"
            printf "ABC neighbors ({wildcards.file}),%s,%s,%s\n" "{params.release}" "{params.base}/{wildcards.file}" "$(date -Iseconds)" > "{log}"
            '''

# --- MERFISH -----------------------------------------------------------------
MER   = ABC.get('merfish', {}) or {}
MER_PRODUCT = MER.get('product')
MER_RELEASE = MER.get('release')
MER_FILES   = MER.get('files') or ({'merfish': MER['file']} if MER.get('file') else {})
MER_FILE_PAT = "(" + "|".join(re.escape(v) for v in MER_FILES.values()) + ")" if (MER_PRODUCT and MER_RELEASE and MER_FILES) else r"^$"
MER_BASE_URL = f"{S3_HOST}/metadata/{MER_PRODUCT}/{MER_RELEASE}" if (MER_PRODUCT and MER_RELEASE) else ""

if MER_PRODUCT and MER_RELEASE and MER_FILES:
    rule abc_merfish_file:
        wildcard_constraints:
            file = MER_FILE_PAT
        params:
            proxy   = WGET_PROXY,
            release = MER_RELEASE,
            base    = MER_BASE_URL
        output:
            out = join(ABC_META, MER_PRODUCT, MER_RELEASE, "{file}")
        log:
            join(ABC_META, MER_PRODUCT, MER_RELEASE, "logs", "{file}.prov.csv")
        shell:
            r'''
            wget {params.proxy} -O "{output.out}" "{params.base}/{wildcards.file}"
            printf "ABC MERFISH ({wildcards.file}),%s,%s,%s\n" "{params.release}" "{params.base}/{wildcards.file}" "$(date -Iseconds)" > "{log}"
            '''

# --- CCF ---------------------------------------------------------------------
CCF   = ABC.get('ccf', {}) or {}
CCF_PRODUCT = CCF.get('product', 'CCF')
CCF_RELEASE = CCF.get('release')
CCF_FILES   = CCF.get('files') or ({'ccf': CCF['file']} if CCF.get('file') else {})
CCF_FILE_PAT = "(" + "|".join(re.escape(v) for v in CCF_FILES.values()) + ")" if (CCF_RELEASE and CCF_FILES) else r"^$"
CCF_BASE_URL = f"{S3_HOST}/metadata/{CCF_PRODUCT}/{CCF_RELEASE}" if CCF_RELEASE else ""

if CCF_RELEASE and CCF_FILES:
    rule abc_ccf_file:
        wildcard_constraints:
            file = CCF_FILE_PAT
        params:
            proxy   = WGET_PROXY,
            release = CCF_RELEASE,
            base    = CCF_BASE_URL
        output:
            out = join(ABC_META, CCF_PRODUCT, CCF_RELEASE, "{file}")
        log:
            join(ABC_META, CCF_PRODUCT, CCF_RELEASE, "logs", "{file}.prov.csv")
        shell:
            r'''
            wget {params.proxy} -O "{output.out}" "{params.base}/{wildcards.file}"
            printf "ABC CCF ({wildcards.file}),%s,%s,%s\n" "{params.release}" "{params.base}/{wildcards.file}" "$(date -Iseconds)" > "{log}"
            '''

def _expand_neighbors():
    return [join(ABC_META, NEIGH_PRODUCT, NEIGH_RELEASE, f)
            for f in NEIGH_FILES.values()] if NEIGH_FILES else []

def _expand_merfish():
    return [join(ABC_META, MER_PRODUCT, MER_RELEASE, f)
            for f in MER_FILES.values()] if (MER_PRODUCT and MER_RELEASE and MER_FILES) else []

def _expand_ccf():
    return [join(ABC_META, CCF_PRODUCT, CCF_RELEASE, f)
            for f in CCF_FILES.values()] if (CCF_RELEASE and CCF_FILES) else []

rule abc_assets:
    input:
        # MapMyCells (mouse & human)
        join(ABC_MAPMY, "mus_musculus", "precomputed_stats.h5"),
        join(ABC_MAPMY, "mus_musculus", "markers.json"),
        join(ABC_MAPMY, "homo_sapiens", "precomputed_stats.h5"),
        join(ABC_MAPMY, "homo_sapiens", "markers.json"),
        # Taxonomy + palettes
        join(ABC_META, "WMB-taxonomy", COLOR_RELEASE, "cluster_annotation_term.csv"),
        ABC_COLORS_LONG,
        ABC_COLORS_JSON,
        # ROI + derived map
        join(ABC_META, ROI_PRODUCT, ROI_RELEASE, "region_of_interest_metadata.csv"),
        join(ABC_DERIVED, "roi_map.tsv")
