# QIAseq region demultiplexing

Microbiome and metagenome use `common/scripts/demultiplex_16s_its.py` through
their existing `rules/filter/scripts/demultiplex_16s_its.py` command-line entry
points. Both Cutadapt rules quote each path and configuration argument. The
shared implementation passes argument lists to Cutadapt and pigz; filesystem
operations and concatenation do not use a shell.

Primer pairs are processed in library-prep configuration order. Each pass takes
only the preceding pass's unmatched pairs, removes its `:region=no_adapter`
marker, and applies the existing paired adapters, overlap/error settings and
20-base minimum length. The second pass removes reverse-complement primers
from short amplicons. Concatenation retains the existing order: final unmatched
pairs first, followed by each configured region. Primer sequences, thresholds
and downstream analysis settings are unchanged.

The CLI uses both `--in1` and `--in2`. Output remains
`<sample>_R1.fastq` and `<sample>_R2.fastq`, or `.fastq.gz` with `--compress`
(pigz level 6, two threads). The sample name is derived from the R1 basename
by removing `_R1.fastq.gz`, as before. The main report is written to
`--log-dir/<sample>_qiaseq_demultiplex.log`, with one appended Cutadapt report
per region. Reverse-complement reports remain under
`--output-dir/rev_comp_log/<sample>_<region>_rev_comp.log`. Reports are replaced
on retry; stderr and nonzero tool exits propagate to the caller. Cutadapt report
contents remain suitable for the existing MultiQC consumer.

Each invocation owns a temporary directory under its output directory. It
ignores old intermediate files and removes only its own temporary directory;
no sample-derived cleanup glob can match another sample. Final FASTQs are
replaced only after both mates have been produced successfully. Each rename is
atomic, but publication of the two files together is not a filesystem
transaction. The workflow must not run the same sample twice concurrently.

Region names must be nonempty, cannot contain path separators, Cutadapt output
template delimiters, `=`, newlines or NUL, and cannot use the reserved name
`unknown`. Empty primer mappings fail before processing. These checks prevent
configuration-derived filenames from escaping the invocation directory.

## Verification

The synthetic regression suite runs real Cutadapt and pigz through both CLI
entry points, plus the actual Snakemake rules without containers. It checks
paired assignment, reverse-complement trimming, unknown markers, empty and
single-region inputs, compressed/uncompressed output, retries, tool failures,
and paths containing spaces, quotes, glob characters and shell metacharacters.

```bash
python -m pytest -q common/.tests/test_qiaseq_demultiplex.py
```

Use an isolated environment with the versions in
`.github/workflows/test-qiaseq.yaml` and pigz on `PATH`. The Python 3.11 tests
use Cutadapt 4.9; operational integration must also exercise the workflow-owned
Cutadapt 3.3 image. Keep the normal microbiome/metagenome `multiqc_report`
dry-runs and the all-workflow server integration gate. Inspect actual QIAseq
logs, region annotations and paired FASTQs, MultiQC reports, delivery archives,
checksums, snapshots and notifications. Synthetic tests do not establish
scientific accuracy or server/container compatibility.
