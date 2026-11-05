"""The R doublet callers are seeded like the Python ones (#297); the scripts are read, not run."""
from pathlib import Path
import re

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / 'rules/quant/scripts'


@pytest.mark.parametrize('script, call', [('scds.R', 'sce_sub <- cxds_bcds_hybrid('), ('scdblfinder.R', 'sce_sub <- scDblFinder(')])
def test_r_doublet_callers_set_the_seed_before_sampling(script, call):
    text = (SCRIPTS / script).read_text()
    assert re.search(r'"--seed",\s*type = "integer",\s*default = 1234', text)
    seeded = text.index('set.seed(args$seed)')
    assert seeded < text.index(call)
    assert 'setup_threads(as.integer(n_threads), args$seed)' in text
    assert text.count('RNGseed = seed') == 2
