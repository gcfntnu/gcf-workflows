"""Configuration-only contract for parent-linked TCR; no server filesystem access."""
from pathlib import PurePosixPath


def validate_singlecell_assay(config):
    """Return the assay, rejecting incomplete TCR settings before WT rule loading."""
    quant = config.get('quant', {})
    assay = quant.get('assay', 'wt')
    if assay not in ('wt', 'tcr'):
        raise ValueError("quant.assay must be 'wt' or 'tcr'")
    if assay == 'wt':
        if config.get('immune'):
            raise ValueError("immune configuration requires quant.assay: tcr")
        return assay
    if quant.get('method') != 'splitpipe':
        raise ValueError('TCR requires quant.method: splitpipe')
    if quant.get('kit') != 'WT' or quant.get('chemistry') != 'v4':
        raise ValueError('This TCR preset supports quant.kit: WT and quant.chemistry: v4')
    geometry = config.get('read_geometry', [])
    if not isinstance(geometry, list) or len(geometry) != 2 or any(type(n) is not int or n <= 0 for n in geometry):
        raise ValueError('TCR requires two positive paired-end read lengths from run metadata')
    if config.get('organism') != 'homo_sapiens':
        raise ValueError('This TCR preset supports organism: homo_sapiens')
    immune = config.get('immune')
    if not isinstance(immune, dict):
        raise ValueError('TCR requires an immune configuration mapping')
    if immune.get('species') != 'human' or immune.get('reference') not in ('default', 'imgt'):
        raise ValueError("TCR requires immune.species: human and reference: default or imgt")
    parent = immune.get('wt_parent')
    if not isinstance(parent, dict):
        raise ValueError('TCR requires immune.wt_parent')
    root = parent.get('root')
    if not isinstance(root, str) or not PurePosixPath(root).is_absolute():
        raise ValueError('immune.wt_parent.root must be an absolute path visible to the execution environment')
    samples = config.get('samples')
    mapping = parent.get('sublibraries')
    if not isinstance(samples, dict) or not samples:
        raise ValueError('TCR requires nonempty config.samples with the TCR sublibrary IDs')
    if not isinstance(mapping, dict) or set(mapping) != set(samples):
        raise ValueError('immune.wt_parent.sublibraries must map every TCR sublibrary exactly once')
    paths = []
    for label, value in list(mapping.items()) + [('combined', parent.get('combined'))]:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'WT parent path for {label} must be a nonempty relative path')
        path = PurePosixPath(value)
        if path.is_absolute() or '..' in path.parts or str(path) == '.':
            raise ValueError(f'WT parent path for {label} must stay beneath immune.wt_parent.root')
        paths.append(str(path))
    if len(set(paths)) != len(paths):
        raise ValueError('WT sublibrary and combined paths must be distinct; parent reuse is not supported')
    order = immune.get('sublibrary_order')
    if not isinstance(order, list) or any(not isinstance(item, str) for item in order):
        raise ValueError('immune.sublibrary_order must be a list of sublibrary IDs')
    if len(order) != len(samples) or set(order) != set(samples):
        raise ValueError('immune.sublibrary_order must contain every TCR sublibrary exactly once')
    if quant.get('use_velo', False) or quant.get('cellbender', {}).get('enabled', False):
        raise ValueError('TCR cannot enable WT velocity or CellBender processing')
    if quant.get('aggregate', {}).get('skip') is not True:
        raise ValueError('TCR requires quant.aggregate.skip: true (WT expression aggregation only)')
    return assay
