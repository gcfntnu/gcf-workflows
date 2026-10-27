"""Resolve WT velocity settings without changing the supplied configuration."""


def resolve_velocity_methods(quant, methods):
    """Apply quantifier defaults, retaining boolean inputs and supporting mixed runs."""
    available = list(methods)
    parse_family = ('parsebio_starsolo', 'parsebio_starsolo_rt')
    if parse_family[0] in available and parse_family[1] not in available:
        available.append(parse_family[1])

    resolved = {method: method != 'splitpipe' for method in available}
    if 'use_velo' not in quant:
        return resolved

    requested = quant['use_velo']
    if type(requested) is bool:
        resolved = {method: requested for method in available}
    elif isinstance(requested, dict):
        for method, enabled in requested.items():
            if method not in available:
                raise ValueError(f"quant.use_velo contains an unselected or unknown method: {method!r}")
            if type(enabled) is not bool:
                raise ValueError(f"quant.use_velo.{method} must be a boolean")
        resolved.update(requested)
        # Both Parse representations share one STARsolo execution configuration.
        parse_values = {requested[method] for method in parse_family if method in requested}
        if len(parse_values) > 1:
            raise ValueError("quant.use_velo settings for parsebio_starsolo and parsebio_starsolo_rt must agree")
        if parse_values:
            enabled = parse_values.pop()
            for method in parse_family:
                if method in resolved:
                    resolved[method] = enabled
    else:
        raise ValueError("quant.use_velo must be a boolean or a mapping of selected quantifiers to booleans")

    if resolved.get('splitpipe', False):
        raise ValueError(
            "Split-Pipe velocity is unsupported. Use parsebio_starsolo for Parse velocity. "
            "Remove quant.use_velo to use quantifier defaults, set it to false, or use a method mapping "
            "with splitpipe: false and parsebio_starsolo: true for a mixed run."
        )
    return resolved
