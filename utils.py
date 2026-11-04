#-*- mode:snakemake -*-
"""setup/configuration, common imports, and  utility functions for all workflows
"""

import warnings
from os.path import join, abspath, dirname, basename
from os import makedirs, environ
import sys
import collections
try:
    from collections import Mapping
except:
    from collections.abc import Mapping
import re
import yaml

try:
    from yaml import CLoader as Loader
except:
    from yaml import Loader as Loader
from datetime import datetime


from snakemake.logging import logger
#from snakemake.workflow import src_gcf
from snakemake.utils import update_config, min_version

min_version("5.10.0")

## setup paths used in workflows
TMPDIR = os.environ.get("TMPDIR", "/tmp")
# config paths will override environment variables
INTERIM_DIR = config.get("interim_dir") or environ.get("GCF_INTERIM", "data/tmp")
makedirs(INTERIM_DIR, exist_ok=True)
EXT_DIR = config.get("ext_dir") or environ.get("GCF_EXT", "data/ext")
makedirs(EXT_DIR, exist_ok=True)
EXT_CACHE = join(EXT_DIR, '.cache')
makedirs(EXT_DIR, exist_ok=True)

FASTQ_DIR = config.get("fastq_dir") or environ.get("GCF_FASTQ", "data/raw/fastq")
while FASTQ_DIR.endswith(os.path.sep):
    FASTQ_DIR = FASTQ_DIR[:-1]
makedirs(FASTQ_DIR, exist_ok=True)

def src_gcf(pth):
    return join(workflow.current_basedir, pth)

GCFDB_DIR = src_gcf("gcfdb")

ORG = config.get('organism', "N/A").lower().strip().replace(' ', '_')
config['organism'] = ORG

PE = len(config['read_geometry']) > 1
BFQ_ALL = []

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
            if isinstance(value, Mapping):
                d[key] = _update(d.get(key, {}), value)
            else:
                if not key in d:
                    try:
                        d[key] = value
                    except:
                        print("Failed dict update in update_config2")
                        print(key)
                        print(d)
                        print(value)
        return d

    return _update(config, extra_config)


def config_val2list(d):
    """Split comma separated values into list of strings.
    """
    for key, value in d.items():
        if isinstance(value, Mapping):
            d[key] = _update(d.get(key, {}), value)
        else:
            for key, val in d.items():
                if isinstance(val, str) and "," in val:
                    d[key] = val.split(",")
        return d
    return _update(config)


default_config_sections = ["db", "quant", "filter", "analysis", "qc", "bfq", "samples", "protocol"]
for section in default_config_sections:
    if section not in config:
        config[section] = {}

# default config
main_fn = src_gcf("main.config")
with open(main_fn) as fh:
    CONF = yaml.load(fh, Loader=Loader) or {}

GCF_SECRET = os.environ.get("GCF_SECRET")
if GCF_SECRET:
    with open(GCF_SECRET) as fh:
        SECRETS = yaml.load(fh, Loader=Loader) or {}
    update_config(CONF, SECRETS)

# library preparation kit specific configuration
libprep_fn = src_gcf("libprep.config")
with open(libprep_fn) as fh:
    LIBPREP_CONF = yaml.load(fh, Loader=Loader) or {}
kit = config.get("libprepkit")
if kit is not None:
    # The installed snakemake command shares bin/ with configmaker.py, which otherwise
    # shadows the configmaker package. Match BFQ's entry-point isolation for this import.
    import_path = sys.path[:]
    executable_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
    try:
        sys.path[:] = [path for path in sys.path if not path or os.path.realpath(path) != executable_dir]
        from configmaker.libprep import LibprepConfig, LibprepConfigError
    finally:
        sys.path[:] = import_path

    # Use the same selection rules as configmaker, including suffix/geometry validation.
    selection = LibprepConfig.load(libprep_fn).select(kit, config['read_geometry'])
    recorded = config.get('libprep_selection')
    if recorded is not None:
        if not isinstance(recorded, Mapping) or recorded.get('entry') != selection.entry:
            raise LibprepConfigError('libprep_selection.entry does not match the selected kit {!r}'.format(selection.entry))
        if recorded.get('read_geometry') != list(selection.read_geometry):
            raise LibprepConfigError('libprep_selection.read_geometry does not match project read_geometry')
        if recorded.get('workflow') != selection.workflow:
            raise LibprepConfigError('libprep_selection.workflow does not match the selected kit workflow')
    kit = selection.entry
    kit_conf = selection.parameters
    # Existing kit-family consumers expect the canonical name without its layout suffix.
    config['libprepkit'] = re.sub(r' (PE|SE)$', '', kit, flags=re.IGNORECASE)
    update_config(CONF, kit_conf)
else:
    if config.get('libprep_selection') is not None:
        raise ValueError('libprep_selection requires libprepkit')
    logger.warning("Running without LIBREPKIT defined!")

# update config (config.yaml). Does not update if key exists
update_config2(config, CONF)


# update read geometry with delta_readlen
if 'delta_readlen' in config and 'read_geometry' in config:
    read_geometry = config["read_geometry"]
    for i, val in enumerate(config['delta_readlen']):
        read_geometry[i] = int(read_geometry[i]) + int(val)
    config["read_geometry"] = read_geometry 


# workflow specific variables 
WORKFLOW = config.get("workflow")
if WORKFLOW is None:
    config["workflow"] = WORKFLOW = "default"
logger.info("WORKFLOW: {}".format(WORKFLOW))
BFQ_INTERIM = join(INTERIM_DIR, WORKFLOW, "bfq")
QC_INTERIM = join(INTERIM_DIR, WORKFLOW, "qc")
QUANT_INTERIM = join(INTERIM_DIR, WORKFLOW, "quant")
FILTER_INTERIM = join(INTERIM_DIR, WORKFLOW, "filter")
ALIGN_INTERIM = join(INTERIM_DIR, WORKFLOW, "align")
ANALYSIS_INTERIM = join(INTERIM_DIR, WORKFLOW, "analysis")


# docker images
docker_fn = src_gcf("docker.config")
with open(docker_fn) as fh:
    DOCKER_CONF = yaml.load(fh, Loader=Loader) or {}
    update_config2(config, DOCKER_CONF)


# load function for statistical models
def load_model(model_yaml_file):
    with open(model_yaml_file) as fh:
        MODELS = yaml.load(fh, Loader=Loader) or {}
        config["models"] = MODELS
        config["model_names"] = list(MODELS.keys())


PATTERNED_FLOWCELL_MACHINES = ('NovaSeq', 'HiSeq 3000', 'HiSeq 4000', 'HiSeq X')
NONPATTERNED_FLOWCELL_MACHINES = ('MiSeq', 'MiniSeq', 'NextSeq', 'HiSeq 2500')

def picard_duplicate_args(machine, remove_optical_duplicates=False):
    """Picard MarkDuplicates options for the sequencing instrument (`machine` from configmaker).

    Duplicates are always tagged (DT:LB or DT:SQ). The optical pixel distance follows Picard's guidance:
    2500 for patterned flowcells, 100 for non-patterned ones and for an unknown (empty) instrument.
    """
    machine = (machine or '').strip()
    if machine.startswith(PATTERNED_FLOWCELL_MACHINES):
        distance = 2500
    elif machine.startswith(NONPATTERNED_FLOWCELL_MACHINES) or not machine:
        distance = 100
    else:
        raise ValueError('machine: no optical-duplicate distance for {!r}; add it to PATTERNED_FLOWCELL_MACHINES '
                         'or NONPATTERNED_FLOWCELL_MACHINES in utils.py'.format(machine))
    args = 'TAGGING_POLICY=All OPTICAL_DUPLICATE_PIXEL_DISTANCE={} '.format(distance)
    if remove_optical_duplicates:
        args += 'REMOVE_SEQUENCING_DUPLICATES=TRUE '
    return args

def get_raw_fastq(wildcards):
    """Returns path to fastq files per sample.
    """
    if 'sublib' in wildcards.keys(): # parsebio: sample==sublib
        R1 = config['samples'][wildcards.sublib].get('R1', [])
        R2 = config['samples'][wildcards.sublib].get('R2', [])
    else:
        R1 = config['samples'][wildcards.sample].get('R1', [])
        R2 = config['samples'][wildcards.sample].get('R2', [])
    if R1 == '':
        R1 = []
    if R2 == '':
        R2 = []
    if R1:
        R1 = R1.split(',')
    if R2:
        R2 = R2.split(',')

    R1 = [join(FASTQ_DIR, i) for i in R1]
    if R2:
        R2 = [join(FASTQ_DIR, i) for i in R2]
        return {'R1': R1, 'R2': R2}
    return {'R1': R1}
