import warnings
from os.path import join, abspath, dirname
from os import makedirs, environ
import sys
import collections
import yaml
try:
    from yaml import CLoader as Loader
except:
    from yaml import Loader as Loader
from datetime import datetime


from snakemake.logging import logger
from snakemake.workflow import srcdir
from snakemake.utils import update_config, min_version

min_version("5.10.0")
TMPDIR = os.environ.get('TMPDIR', '/tmp')
# environment variables can override config file
INTERIM_DIR = environ.get('GCF_INTERIM') or config.get('interim_dir', 'data/tmp/')
makedirs(INTERIM_DIR, exist_ok=True)
EXT_DIR = environ.get('GCF_EXT') or config.get('ext_dir', 'data/ext')
FASTQ_DIR =  environ.get('GCF_FASTQ') or config.get('fastq_dir','data/raw/fastq')
makedirs(FASTQ_DIR, exist_ok=True)

if not 'db' in config:
    config['db'] = {}

def update_config2(config, extra_config):
    """Recursively update dictionary config with overwrite_config.

    From Snakemake codebase (update_config), this does not overwrite if key exist

    See
    http://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    for details.

    Args:
      config (dict): dictionary to update
      overwrite_config (dict): dictionary whose items will overwrite those in config

    """

    def _update(d, u):
        for (key, value) in u.items():
            if (isinstance(value, collections.Mapping)):
                d[key] = _update(d.get(key, {}), value)
            else:
                if not key in d:
                    d[key] = value
        return d
    return _update(config, extra_config)


default_config_sections = ['db', 'quant', 'filter', 'analysis', 'samples']
for section in default_config_sections:
    if section not in config:
        config[section] = {}

# load function for statistical models
def load_model(model_yaml_file):
    with open(model_yaml_file) as fh:
        MODELS  = yaml.load(fh, Loader=Loader) or {}
        config['models'] = MODELS
        config['model_names'] = list(MODELS.keys())

# library preparation kit specific configuration
libprep_fn = srcdir('libprep.config')
with open(libprep_fn) as fh:
    LIBPREP_CONF  = yaml.load(fh, Loader=Loader) or {}
kit = config.get('libprepkit')
if kit in LIBPREP_CONF:
    LIBPREP = LIBPREP_CONF[kit]
    if 'reference_db' in LIBPREP:
        config['db']['reference_db'] = LIBPREP['reference_db']
else:
    if kit is None:
        logger.warning('Running without LIBREPKIT defined!')
    else:
        logger.warning('`{}` is not a valid librepkit name'.format(kit))
        sys.exit()
    

# docker images
docker_fn = srcdir('docker.config')
with open(docker_fn) as fh:
    dck = yaml.load(fh, Loader=Loader) or {}
    update_config2(config, dck)
