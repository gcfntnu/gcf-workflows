"""Load single-cell defaults and diagnose required settings without changing overrides."""

from collections.abc import Mapping

import yaml


def singlecell_config_error(config, defaults_path, message):
    quant = config.get("quant", {})
    method = quant.get("method", "<not selected>") if isinstance(quant, Mapping) else "<not selected>"
    return ValueError(f"Single-cell configuration for quant.method={method!r}: {message}. Defaults file: {defaults_path}")


def load_singlecell_defaults(config, defaults_path):
    try:
        with open(defaults_path) as handle:
            defaults = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as error:
        raise singlecell_config_error(config, defaults_path, f"cannot load required defaults ({error})") from error
    if not isinstance(defaults, Mapping) or not defaults:
        raise singlecell_config_error(config, defaults_path, "required defaults must be a nonempty YAML mapping")
    return defaults


def required_singlecell_setting(config, key, defaults_path):
    value = config
    for part in key.split("."):
        if not isinstance(value, Mapping) or part not in value or value[part] is None:
            raise singlecell_config_error(config, defaults_path, f"missing required setting {key!r}")
        value = value[part]
    return value
