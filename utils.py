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
GCFDB_DIR = srcdir("gcfdb")


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
            if isinstance(value, collections.Mapping):
                d[key] = _update(d.get(key, {}), value)
            else:
                if not key in d:
                    try:
                        d[key] = value
                    except:
                        print("MMMM")
                        print(key)
                        print(d)
        return d

    return _update(config, extra_config)


def config_val2list(d):
    """Split comma separated values into list of strings.
    """
    for key, value in d.items():
        if isinstance(value, collections.Mapping):
            d[key] = _update(d.get(key, {}), value)
        else:
            for key, val in d.items():
                if isinstance(val, str) and "," in val:
                    d[key] = val.split(",")
        return d
    return _update(config)


default_config_sections = ["db", "quant", "filter", "analysis", "qc", "bfq", "samples"]
for section in default_config_sections:
    if section not in config:
        config[section] = {}

# default config
main_fn = srcdir("main.config")
with open(main_fn) as fh:
    CONF = yaml.load(fh, Loader=Loader) or {}

GCF_SECRET = os.environ.get("GCF_SECRET")
if GCF_SECRET:
    with open(GCF_SECRET) as fh:
        SECRETS = yaml.load(fh, Loader=Loader) or {}
    update_config(CONF, SECRETS)

# library preparation kit specific configuration
libprep_fn = srcdir("libprep.config")
with open(libprep_fn) as fh:
    LIBPREP_CONF = yaml.load(fh, Loader=Loader) or {}
kit = config.get("libprepkit")
if kit is not None:
    if len(config["read_geometry"]) > 1:
        kit += " PE"
    else:
        kit += " SE"
if kit in LIBPREP_CONF:
    # overwrite default config
    update_config(CONF, LIBPREP_CONF[kit])
else:
    if kit is None:
        logger.warning("Running without LIBREPKIT defined!")
    else:
        logger.warning("`{}` is not a valid librepkit name".format(kit))
        sys.exit()

# update config (config.yaml). Does not update if key exists
update_config2(config, CONF)

# workflow specific variables 
WORKFLOW = config.get("workflow")
if WORKFLOW is None:
    config["workflow"] = WORKFLOW = "default"
BFQ_INTERIM = join(INTERIM_DIR, WORKFLOW, "bfq")
QC_INTERIM = join(INTERIM_DIR, WORKFLOW, "qc")
QUANT_INTERIM = join(INTERIM_DIR, WORKFLOW, "quant")
FILTER_INTERIM = join(INTERIM_DIR, WORKFLOW, "filter")
ALIGN_INTERIM = join(INTERIM_DIR, WORKFLOW, "align")
ANALYSIS_INTERIM = join(INTERIM_DIR, WORKFLOW, "analysis")



# docker images
docker_fn = srcdir("docker.config")
with open(docker_fn) as fh:
    DOCKER_CONF = yaml.load(fh, Loader=Loader) or {}
    update_config2(config, DOCKER_CONF)


# load function for statistical models
def load_model(model_yaml_file):
    with open(model_yaml_file) as fh:
        MODELS = yaml.load(fh, Loader=Loader) or {}
        config["models"] = MODELS
        config["model_names"] = list(MODELS.keys())


def get_raw_fastq(wildcards):
    """Returns path to fastq files per sample.
    """

    R1 = config["samples"][wildcards.sample].get("R1", [])
    R2 = config["samples"][wildcards.sample].get("R2", [])
    if R1 == "":
        R1 = []
    if R2 == "":
        R2 = []
    if R1:
        R1 = R1.split(",")
    if R2:
        R2 = R2.split(",")

    R1 = [join(FASTQ_DIR, i) for i in R1]
    if R2:
        R2 = [join(FASTQ_DIR, i) for i in R2]
        return {"R1": R1, "R2": R2}
    else:
        return {"R1": R1}


def get_processed_fastq(wildcards):
    ZIP_FILTERED_FASTQ = config["filter"].get("zip_filtered_fastq", False)
    DST_PTH = join(FILTER_INTERIM, "cleaned")
    fastq = get_raw_fastq(wildcards)
    R1 = [i.split(FASTQ_DIR)[-1][1:] for i in fastq["R1"]]
    if not ZIP_FILTERED_FASTQ:
        R1 = [i.split(".gz")[0] for i in R1 if i.endswith(".gz")]
    R2 = fastq.get("R2")
    if R2:
        R2 = [i.split(FASTQ_DIR)[-1][1:] for i in fastq["R2"]]
        if not ZIP_FILTERED_FASTQ:
            R2 = [i.split(".gz")[0] for i in R2 if i.endswith(".gz")]
    if ZIP_FILTERED_FASTQ:
        ext = ".fastq.gz"
    else:
        ext = ".fastq"
    R1 = [join(DST_PTH, wildcards.sample + "_R1" + ext)]
    if R2:
        R2 = [join(DST_PTH, wildcards.sample + "_R2" + ext)]
        out = {"R1": R1, "R2": R2}
    else:
        out = {"R1": R1}
    return out


# single cell get fastq functions
def umitools_get_filtered_fastq_R1(wildcards):
    if config["filter"]["trim"]["quantifier"] == "skip":
        return get_raw_fastq(wildcards)["R1"]
    return [join(FILTER_INTERIM, "cleaned", "{}_R1.fastq.gz".format(wildcards.sample))]


def umitools_get_filtered_fastq_R2(wildcards):
    if config["filter"]["trim"]["quantifier"] == "skip":
        return get_raw_fastq(wildcards)["R2"]
    return [join(FILTER_INTERIM, "cleaned", "{}_R2.fastq.gz".format(wildcards.sample))]


def get_filtered_fastq_R1(wildcards):
    if config["filter"]["trim"]["quantifier"] == "skip":
        return get_raw_fastq(wildcards)["R1"]
    elif config["libprepkit"].startswith("Drop-Seq"):
        return umitools_get_filtered_fastq_R1(wildcards)
    return [
        join(
            FILTER_INTERIM,
            "cleaned",
            "{}_S1_L000_R1_001.fastq.gz".format(wildcards.sample),
        )
    ]


def get_filtered_fastq_R2(wildcards):
    if config["filter"]["trim"]["quantifier"] == "skip":
        return get_raw_fastq(wildcards)["R2"]
    elif config["libprepkit"].startswith("Drop-Seq"):
        return umitools_get_filtered_fastq_R2(wildcards)
    return [
        join(
            FILTER_INTERIM,
            "cleaned",
            "{}_S1_L000_R2_001.fastq.gz".format(wildcards.sample),
        )
    ]
