# Velocity settings for single-cell workflows

WT velocity settings are resolved for each selected quantifier. Standard Split-Pipe
expression processing runs without velocity generation or attachment. Use Parse
STARsolo when Parse velocity matrices are required.

| Quantifier | Default when `quant.use_velo` is omitted | Explicit velocity request |
| --- | --- | --- |
| `splitpipe` | Disabled | Rejected before workflow execution |
| `parsebio_starsolo` / `parsebio_starsolo_rt` | Enabled | Existing STARsolo matrix generation and attachment |
| `10x_starsolo` | Enabled | Existing STARsolo matrix generation and attachment |

Other quantifiers retain their previous default flag behavior; this change adds no
velocity capability to them. Parent-linked Parse TCR retains its existing explicit
`use_velo: false` setting and assay validation.

## Configuration

An ordinary Split-Pipe WT configuration needs no velocity override:

```yaml
quant:
  method: splitpipe
```

The workflow applies method defaults after loading the project, library-prep and
shared configuration. The shared defaults no longer insert `use_velo: true` into
every workflow. An omitted setting resolves to false for Split-Pipe and preserves
STARsolo's true default, including when the quantifier is changed from Split-Pipe
to Parse STARsolo in a project configuration.

Existing boolean settings remain accepted. `use_velo: false` disables velocity for
all selected methods. `use_velo: true` requests it for all selected methods and is
therefore rejected if Split-Pipe is selected. The error identifies the unsupported
route and directs the operator to Parse STARsolo.

WT configurations can also supply a mapping for independent method selection:

```yaml
quant:
  method: splitpipe,parsebio_starsolo
  use_velo:
    splitpipe: false
    parsebio_starsolo: true
```

This keeps Split-Pipe expression processing free of velocity work while retaining
Parse STARsolo velocity generation, dependencies and AnnData layers. Omitting the
mapping produces the same defaults. A mapping entry takes a YAML boolean; omitted
entries use the method default. Unknown or unselected method names and nonboolean
values fail during workflow loading.

The automatically available `parsebio_starsolo_rt` representation shares the
`parsebio_starsolo` velocity choice. Setting either representation also sets the
other when both are available; explicit conflicting values are rejected. This
keeps their shared STARsolo execution configuration consistent.

The existing BFQ delivery-method selection is unchanged. Selecting multiple
quantifiers does not promise that every quantifier's results are included in a
standard BFQ bundle; request the intended analysis targets as usual.

## Existing configurations and retained workdirs

Normal WT configurations produced by configmaker omit the old shared velocity
default, so they acquire the corrected behavior when loaded by this workflow.
Existing explicit `use_velo: false` configurations continue to work.

An old saved or manually expanded configuration may contain `use_velo: true`.
There is no reliable provenance distinguishing a copied historical default from a
deliberate request. If it selects Split-Pipe, review the intent and either remove
the setting to use method defaults, set it to false, or use the mapping above.
The workflow reports the unsupported request instead of silently discarding it.

Retained BFQ analysis workdirs contain their own workflow and configuration copies.
Updating the installed workflow alone does not change those copies. Apply reviewed
configuration/source changes through the normal operator recovery procedure;
this change does not rewrite retained state, inputs or outputs.

## Support boundary and verification

The custom `splitpipe_velo.py` helper derives spliced/unspliced matrices from
transcript-assignment tables. It remains available, but is not a supported BFQ
biological-sample velocity delivery route. The matrix converter rejects explicit
Split-Pipe velocity attachment before reading expression matrices. Full Split-Pipe
velocity generation/aggregation and BAM-based reconstruction are outside this change.

`singlecell/.tests/test_velocity.py` exercises actual workflow configuration merging,
WT presets, STARsolo outputs/arguments, mixed-method aggregation, unsupported
requests and direct converter validation using synthetic inputs and dry runs.
Run it with the dependencies in `singlecell/.tests/requirements.txt`; it executes
no quantifiers, containers, downloads or SMTP delivery.

Before merge, include a normal Parse Split-Pipe run without a velocity override in
the arranged server integration campaign. Check expression matrices, sample/library
metadata, doublet detection, QC, preprocessing, reports and standard delivery.
Where STARsolo velocity is exercised, check its declared matrices and attached
layers using the existing scientific validation procedure. Local dry runs establish
configuration and DAG behavior, not scientific output validity.
