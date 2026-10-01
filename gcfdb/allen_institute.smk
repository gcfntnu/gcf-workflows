# allen_institute.smk — Allen Brain Cell Atlas downloads and derived resources
import re
from os.path import join
from urllib.parse import quote as _q


# -----------------------------------------------------------------------------
# Config and output roots
# -----------------------------------------------------------------------------
ABC = config['db']['allen_abc']
S3 = ABC.get('s3', {})
BUCKET = S3.get('bucket', 'allen-brain-cell-atlas').strip('/')
REGION = S3.get('region', 'us-west-2')
S3_HOST = f"https://{BUCKET}.s3.{REGION}.amazonaws.com"

WGET_PROXY = config.get('proxy', {}).get('wget', '')
DOCKER_DEF = 'docker://' + config['docker']['default']

EXT_ABC = join(EXT_DIR, 'allen-brain-cell-atlas')
ABC_MAPMY = join(EXT_ABC, 'mapmycells')
ABC_META = join(EXT_ABC, 'metadata')
ABC_DERIVED = join(EXT_ABC, 'derived')

MAPMYCELLS = ABC['mapmycells']
TAXONOMY = ABC['taxonomy']
TAXONOMY_MOUSE_ADDON = ABC.get('taxonomy_mouse_addon', {})

CANONICAL_TAXONOMY_FILES = {
    'cluster': 'cluster.csv',
    'term': 'cluster_annotation_term.csv',
    'term_set': 'cluster_annotation_term_set.csv',
    'membership': 'cluster_to_cluster_annotation_membership.csv',
}

SUPPORTED_TAXONOMY_ORGANISMS = tuple(TAXONOMY.keys())
TAXONOMY_PRODUCTS = tuple(TAXONOMY[organism]['product'] for organism in SUPPORTED_TAXONOMY_ORGANISMS)
TAXONOMY_RELEASES = tuple(TAXONOMY[organism]['release'] for organism in SUPPORTED_TAXONOMY_ORGANISMS)


# -----------------------------------------------------------------------------
# URL and path helpers
# -----------------------------------------------------------------------------
def url_mapmycells(organism, kind):
    """Return mapmycells/<product>/<release>/<file> URL for an organism."""
    spec = MAPMYCELLS[organism]
    key = f"mapmycells/{spec['product']}/{spec['release']}/{spec['files'][kind]}"
    return f"{S3_HOST}/{_q(key, safe='/')}"


def url_metadata(product, release, filename):
    """Return metadata/<product>/<release>/<file> URL."""
    key = f"metadata/{product}/{release}/{filename}"
    return f"{S3_HOST}/{_q(key, safe='/')}"


def abc_taxonomy_file(organism, kind):
    """Return the local path to a canonical taxonomy table."""
    if organism not in TAXONOMY:
        raise ValueError(f"No Allen taxonomy configured for organism: {organism}")
    if kind not in CANONICAL_TAXONOMY_FILES:
        raise ValueError(f"Unknown canonical Allen taxonomy file kind: {kind}")

    spec = TAXONOMY[organism]
    return join(
        ABC_META,
        spec['product'],
        spec['release'],
        CANONICAL_TAXONOMY_FILES[kind],
    )


def abc_taxonomy_colors(organism, fmt='json'):
    """Return the derived Allen taxonomy color palette for an organism."""
    if organism not in TAXONOMY:
        raise ValueError(f"No Allen taxonomy configured for organism: {organism}")
    if fmt not in ('json', 'long'):
        raise ValueError(f"Unknown Allen taxonomy color format: {fmt}")

    filename = 'abc_colors.json' if fmt == 'json' else 'abc_colors.long.tsv'
    return join(ABC_DERIVED, 'taxonomy', organism, filename)


def abc_mouse_taxonomy_addon_file(kind):
    """Return a configured mouse-only taxonomy addon file."""
    files = TAXONOMY_MOUSE_ADDON.get('files', {})
    if kind not in files:
        raise ValueError(f"Unknown Allen mouse taxonomy addon file kind: {kind}")

    return join(
        ABC_META,
        TAXONOMY_MOUSE_ADDON['product'],
        TAXONOMY_MOUSE_ADDON['release'],
        files[kind],
    )


def _taxonomy_spec_from_product(product):
    matches = [spec for spec in TAXONOMY.values() if spec['product'] == product]
    if len(matches) != 1:
        raise ValueError(f"Allen taxonomy product is not uniquely configured: {product}")
    return matches[0]


def _canonical_taxonomy_url(wildcards):
    spec = _taxonomy_spec_from_product(wildcards.taxonomy_product)
    if wildcards.taxonomy_release != spec['release']:
        raise ValueError(
            f"Unexpected release for {wildcards.taxonomy_product}: "
            f"{wildcards.taxonomy_release} != {spec['release']}"
        )
    if wildcards.taxonomy_file not in CANONICAL_TAXONOMY_FILES.values():
        raise ValueError(f"Not a canonical Allen taxonomy file: {wildcards.taxonomy_file}")

    return url_metadata(
        wildcards.taxonomy_product,
        wildcards.taxonomy_release,
        wildcards.taxonomy_file,
    )


def _taxonomy_levels(wildcards):
    return ' '.join(TAXONOMY[wildcards.organism]['levels'])


# -----------------------------------------------------------------------------
# MapMyCells reference assets
# -----------------------------------------------------------------------------
rule mapmycells_precomputed_mouse:
    params:
        url = url_mapmycells('mus_musculus', 'precomputed'),
        release = MAPMYCELLS['mus_musculus']['release'],
        name = 'ABC precomputed (mouse)',
        proxy = WGET_PROXY
    output:
        pre_stats = join(ABC_MAPMY, 'mus_musculus', 'precomputed_stats.h5')
    log:
        join(ABC_MAPMY, 'mus_musculus', 'logs', 'precomputed.prov.csv')
    shell:
        r'''
        wget {params.proxy} -O "{output.pre_stats}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''


rule mapmycells_markers_mouse:
    params:
        url = url_mapmycells('mus_musculus', 'markers'),
        release = MAPMYCELLS['mus_musculus']['release'],
        name = 'ABC markers (mouse)',
        proxy = WGET_PROXY
    output:
        markers = join(ABC_MAPMY, 'mus_musculus', 'markers.json')
    log:
        join(ABC_MAPMY, 'mus_musculus', 'logs', 'markers.prov.csv')
    shell:
        r'''
        wget {params.proxy} -O "{output.markers}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''


rule mapmycells_precomputed_human:
    params:
        url = url_mapmycells('homo_sapiens', 'precomputed'),
        release = MAPMYCELLS['homo_sapiens']['release'],
        name = 'ABC precomputed (human)',
        proxy = WGET_PROXY
    output:
        pre_stats = join(ABC_MAPMY, 'homo_sapiens', 'precomputed_stats.h5')
    log:
        join(ABC_MAPMY, 'homo_sapiens', 'logs', 'precomputed.prov.csv')
    shell:
        r'''
        wget {params.proxy} -O "{output.pre_stats}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''


rule mapmycells_markers_human:
    params:
        url = url_mapmycells('homo_sapiens', 'markers'),
        release = MAPMYCELLS['homo_sapiens']['release'],
        name = 'ABC markers (human)',
        proxy = WGET_PROXY
    output:
        markers = join(ABC_MAPMY, 'homo_sapiens', 'markers.json')
    log:
        join(ABC_MAPMY, 'homo_sapiens', 'logs', 'markers.prov.csv')
    shell:
        r'''
        wget {params.proxy} -O "{output.markers}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''


# -----------------------------------------------------------------------------
# Canonical taxonomy assets: WMB and WHB
# -----------------------------------------------------------------------------
rule abc_taxonomy_file:
    wildcard_constraints:
        taxonomy_product = '|'.join(re.escape(product) for product in TAXONOMY_PRODUCTS),
        taxonomy_release = '|'.join(re.escape(release) for release in TAXONOMY_RELEASES),
        taxonomy_file = '|'.join(re.escape(filename) for filename in CANONICAL_TAXONOMY_FILES.values())
    params:
        url = _canonical_taxonomy_url,
        proxy = WGET_PROXY
    output:
        out = join(ABC_META, '{taxonomy_product}', '{taxonomy_release}', '{taxonomy_file}')
    log:
        join(ABC_META, '{taxonomy_product}', '{taxonomy_release}', 'logs', '{taxonomy_file}.prov.csv')
    shell:
        r'''
        wget {params.proxy} -O "{output.out}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{wildcards.taxonomy_product}" "{wildcards.taxonomy_release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''


rule abc_build_colors:
    wildcard_constraints:
        organism = '|'.join(re.escape(organism) for organism in SUPPORTED_TAXONOMY_ORGANISMS)
    input:
        term_csv = lambda wc: abc_taxonomy_file(wc.organism, 'term')
    output:
        long = join(ABC_DERIVED, 'taxonomy', '{organism}', 'abc_colors.long.tsv'),
        json = join(ABC_DERIVED, 'taxonomy', '{organism}', 'abc_colors.json')
    params:
        script = src_gcf('gcfdb/scripts/build_abc_colors_from_taxonomy.py'),
        levels = _taxonomy_levels
    container:
        DOCKER_DEF
    shell:
        'python {params.script} '
        '--term-csv {input.term_csv} '
        '--out-long {output.long} '
        '--out-json {output.json} '
        '--levels {params.levels} '


# -----------------------------------------------------------------------------
# Mouse-only taxonomy additions
# -----------------------------------------------------------------------------
if TAXONOMY_MOUSE_ADDON:
    rule abc_mouse_taxonomy_cluster_metadata:
        params:
            url = url_metadata(
                TAXONOMY_MOUSE_ADDON['product'],
                TAXONOMY_MOUSE_ADDON['release'],
                TAXONOMY_MOUSE_ADDON['files']['cluster_metadata'],
            ),
            release = TAXONOMY_MOUSE_ADDON['release'],
            name = 'ABC mouse taxonomy cluster metadata',
            proxy = WGET_PROXY
        output:
            metadata_xlsx = abc_mouse_taxonomy_addon_file('cluster_metadata')
        log:
            join(
                ABC_META,
                TAXONOMY_MOUSE_ADDON['product'],
                TAXONOMY_MOUSE_ADDON['release'],
                'logs',
                'cluster_metadata.prov.csv',
            )
        shell:
            r'''
            wget {params.proxy} -O "{output.metadata_xlsx}" "{params.url}"
            printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
            '''


# -----------------------------------------------------------------------------
# Mouse atlas spatial/anatomical resources
# -----------------------------------------------------------------------------
ROI = ABC.get('roi', {})
ROI_PRODUCT = ROI.get('product', MAPMYCELLS['mus_musculus']['product'])
ROI_RELEASE = ROI.get('release', MAPMYCELLS['mus_musculus']['release'])

rule abc_roi_metadata:
    params:
        url = url_metadata(ROI_PRODUCT, ROI_RELEASE, 'region_of_interest_metadata.csv'),
        release = ROI_RELEASE,
        name = f'ABC ROI metadata ({ROI_PRODUCT})',
        proxy = WGET_PROXY
    output:
        roi_csv = join(ABC_META, ROI_PRODUCT, ROI_RELEASE, 'region_of_interest_metadata.csv')
    log:
        join(ABC_META, ROI_PRODUCT, ROI_RELEASE, 'logs', 'roi_metadata.prov.csv')
    shell:
        r'''
        wget {params.proxy} -O "{output.roi_csv}" "{params.url}"
        printf "%s,%s,%s,%s\n" "{params.name}" "{params.release}" "{params.url}" "$(date -Iseconds)" > "{log}"
        '''


rule abc_build_roi_map:
    input:
        roi = rules.abc_roi_metadata.output.roi_csv
    output:
        out = join(ABC_DERIVED, 'roi_map.tsv')
    params:
        script = src_gcf('scripts/build_abc_roi_map.py')
    container:
        DOCKER_DEF
    shell:
        'python {params.script} '
        '--roi {input.roi} '
        '--out {output.out} '


NEIGH = ABC.get('neighbors', {}) or {}
NEIGH_PRODUCT = NEIGH.get('product', MAPMYCELLS['mus_musculus']['product'])
NEIGH_RELEASE = NEIGH.get('release', ROI_RELEASE)
NEIGH_FILES = NEIGH.get('files') or ({'neighbors': NEIGH['file']} if NEIGH.get('file') else {})
NEIGH_FILE_PAT = '|'.join(re.escape(v) for v in NEIGH_FILES.values()) if NEIGH_FILES else r'^$'
NEIGH_BASE_URL = f"{S3_HOST}/metadata/{NEIGH_PRODUCT}/{NEIGH_RELEASE}"

if NEIGH_FILES:
    rule abc_neighbors_file:
        wildcard_constraints:
            file = NEIGH_FILE_PAT
        params:
            proxy = WGET_PROXY,
            release = NEIGH_RELEASE,
            base = NEIGH_BASE_URL
        output:
            out = join(ABC_META, NEIGH_PRODUCT, NEIGH_RELEASE, '{file}')
        log:
            join(ABC_META, NEIGH_PRODUCT, NEIGH_RELEASE, 'logs', '{file}.prov.csv')
        shell:
            r'''
            wget {params.proxy} -O "{output.out}" "{params.base}/{wildcards.file}"
            printf "ABC neighbors ({wildcards.file}),%s,%s,%s\n" "{params.release}" "{params.base}/{wildcards.file}" "$(date -Iseconds)" > "{log}"
            '''


MER = ABC.get('merfish', {}) or {}
MER_PRODUCT = MER.get('product')
MER_RELEASE = MER.get('release')
MER_FILES = MER.get('files') or ({'merfish': MER['file']} if MER.get('file') else {})
MER_FILE_PAT = '|'.join(re.escape(v) for v in MER_FILES.values()) if MER_FILES else r'^$'
MER_BASE_URL = f"{S3_HOST}/metadata/{MER_PRODUCT}/{MER_RELEASE}" if (MER_PRODUCT and MER_RELEASE) else ''

if MER_PRODUCT and MER_RELEASE and MER_FILES:
    rule abc_merfish_file:
        wildcard_constraints:
            file = MER_FILE_PAT
        params:
            proxy = WGET_PROXY,
            release = MER_RELEASE,
            base = MER_BASE_URL
        output:
            out = join(ABC_META, MER_PRODUCT, MER_RELEASE, '{file}')
        log:
            join(ABC_META, MER_PRODUCT, MER_RELEASE, 'logs', '{file}.prov.csv')
        shell:
            r'''
            wget {params.proxy} -O "{output.out}" "{params.base}/{wildcards.file}"
            printf "ABC MERFISH ({wildcards.file}),%s,%s,%s\n" "{params.release}" "{params.base}/{wildcards.file}" "$(date -Iseconds)" > "{log}"
            '''


CCF = ABC.get('ccf', {}) or {}
CCF_PRODUCT = CCF.get('product', 'CCF')
CCF_RELEASE = CCF.get('release')
CCF_FILES = CCF.get('files') or ({'ccf': CCF['file']} if CCF.get('file') else {})
CCF_FILE_PAT = '|'.join(re.escape(v) for v in CCF_FILES.values()) if CCF_FILES else r'^$'
CCF_BASE_URL = f"{S3_HOST}/metadata/{CCF_PRODUCT}/{CCF_RELEASE}" if CCF_RELEASE else ''

if CCF_RELEASE and CCF_FILES:
    rule abc_ccf_file:
        wildcard_constraints:
            file = CCF_FILE_PAT
        params:
            proxy = WGET_PROXY,
            release = CCF_RELEASE,
            base = CCF_BASE_URL
        output:
            out = join(ABC_META, CCF_PRODUCT, CCF_RELEASE, '{file}')
        log:
            join(ABC_META, CCF_PRODUCT, CCF_RELEASE, 'logs', '{file}.prov.csv')
        shell:
            r'''
            wget {params.proxy} -O "{output.out}" "{params.base}/{wildcards.file}"
            printf "ABC CCF ({wildcards.file}),%s,%s,%s\n" "{params.release}" "{params.base}/{wildcards.file}" "$(date -Iseconds)" > "{log}"
            '''


def _expand_neighbors():
    return [join(ABC_META, NEIGH_PRODUCT, NEIGH_RELEASE, f) for f in NEIGH_FILES.values()] if NEIGH_FILES else []


def _expand_merfish():
    return [join(ABC_META, MER_PRODUCT, MER_RELEASE, f) for f in MER_FILES.values()] if MER_FILES else []


def _expand_ccf():
    return [join(ABC_META, CCF_PRODUCT, CCF_RELEASE, f) for f in CCF_FILES.values()] if CCF_FILES else []


def _canonical_taxonomy_assets():
    return [
        abc_taxonomy_file(organism, kind)
        for organism in SUPPORTED_TAXONOMY_ORGANISMS
        for kind in CANONICAL_TAXONOMY_FILES
    ]


def _taxonomy_color_assets():
    return [
        path
        for organism in SUPPORTED_TAXONOMY_ORGANISMS
        for path in (abc_taxonomy_colors(organism, 'long'), abc_taxonomy_colors(organism, 'json'))
    ]


def _mouse_taxonomy_addons():
    if not TAXONOMY_MOUSE_ADDON:
        return []
    return [abc_mouse_taxonomy_addon_file(kind) for kind in TAXONOMY_MOUSE_ADDON.get('files', {})]


rule abc_assets:
    input:
        join(ABC_MAPMY, 'mus_musculus', 'precomputed_stats.h5'),
        join(ABC_MAPMY, 'mus_musculus', 'markers.json'),
        join(ABC_MAPMY, 'homo_sapiens', 'precomputed_stats.h5'),
        join(ABC_MAPMY, 'homo_sapiens', 'markers.json'),
        _canonical_taxonomy_assets(),
        _taxonomy_color_assets(),
        _mouse_taxonomy_addons(),
        join(ABC_META, ROI_PRODUCT, ROI_RELEASE, 'region_of_interest_metadata.csv'),
        join(ABC_DERIVED, 'roi_map.tsv'),
        _expand_neighbors(),
        _expand_merfish(),
        _expand_ccf()
