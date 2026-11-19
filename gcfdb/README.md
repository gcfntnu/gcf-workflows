# gcfdb
Snakemake rules for downloading data from external sources

## Resource identity and test overrides

Resource rules resolve paths below `EXT_DIR`, selected by `config.ext_dir`,
`GCF_EXT`, or the default `data/ext` (see [utils.py](../utils.py)). A test setup
can reuse an existing resource tree while selecting a distinct release or
assembly. Keep reduced resources beside full resources under a different
identity; do not replace a full reference or database with test contents.

The values below describe configuration interfaces. Actual site selections,
private fixtures, provenance and storage locations belong in the operator's
private configuration/runbook. Changing a test override does not change the
suite's scientific defaults.

### Project reference for Cell Ranger

For a genome augmented by a reporter sequence, point the project's `config.yaml` at the prepared directory.

```yaml
ref_dir: refs/panel-v1
```

The directory must contain the augmented FASTA and its matching annotation. Relative paths start at the project working directory.

| File | Contents |
| --- | --- |
| `refs/panel-v1/fasta/genome.fa` | Genome sequences including the reporter contig |
| `refs/panel-v1/anno/genes.gtf` | Gene annotation including the reporter feature with the same contig name and coordinates |

The builder copies the selected FASTA and GTF into a private working directory, verifies their checksums and runs
`cellranger mkref --fasta refs/panel-v1/fasta/genome.fa --genes refs/panel-v1/anno/genes.gtf` there.
Absolute project paths and paths containing `..` are mapped into that private directory. The published index retains the verified originals in `gcf-source/`.
Count passes its workdir’s `index` symlink to `--transcriptome`. The link pins the published generation in `refs/panel-v1/cellranger/builds/<source-digest>/<build-key>`.
When `ref_dir` is a symlink, the workflow resolves the selected directory before building. Changing the symlink later cannot redirect a planned count.
The source digest covers the SHA-256 checksums of both files. The build key covers the normalized Cell Ranger genome identifier, configured image and recipe version.
Workdir receipts distinguish the canonical cache directory, so changing to another prepared directory selects its own index even if the inputs are identical.
Each complete index is published with a manifest recording every file's size and SHA-256 checksum.
When a target uses the pinned generation, planning in a retained workdir checks its files against that manifest, including on dry runs.
This reads the full published index and its frozen FASTA and GTF once per check, which can add substantial I/O for a large STAR index.
A forced publisher or a second project verifies the generation before pinning it for count. A workdir pin can be reused only for its original generation.
If a managed source is missing, the workflow fetches it and creates a fresh local pin. Once both sources exist, a rebuild must match the input digest from planning.
Changing either source, the genome identifier or the image selects a separate directory. Existing counts can continue reading their original index.
Legacy indexes directly under `ref_dir/cellranger` lack this identity and require a fresh build.
The report derives its gene metadata and exon sequences in the project workdir under `<QUANT_INTERIM>/reference/<source-digest>/` from the same input pair.
The directory name is arbitrary. The configured provider and annotation settings still apply.

With the standard 10x kit and no project reference directory, the workflow selects the reference by organism.

- For human, mouse and the combined human–mouse reference, it uses the supported 10x prebuilt.
  Omit `db.10xgenomics.release` in the project to inherit the configured default. An empty or null effective release stops planning.
- For an organism with an Ensembl reference and no supported prebuilt, the workflow selects Ensembl for the build, count and gene metadata.
  For example, `salmo_salar` at Ensembl release 116 uses assembly `Ssal_v3.1` under `<EXT_DIR>/ensembl/release-116/salmo_salar/Ssal_v3.1/`.
  For organisms in an Ensembl Genomes division, the workflow uses `db.ensemblgenomes.release` for the archives and reports the division.
  For example, `arabidopsis_thaliana` uses the release 56 Ensembl Plants FASTA and annotation under
  `<EXT_DIR>/ensemblgenomes/release-56/plants/arabidopsis_thaliana/TAIR10/`.
  Projects selecting `db.reference_db: ensembl` and projects using the 10x fallback share this division- and release-specific cache.
  On first use, the workflow fetches these archives even if the older `ensembl/release-<main-release>/` layout contains a nonvertebrate genome. That path cannot establish the Ensembl Genomes release. An explicit `ref_dir` can still select a prepared directory there.
  These routes and the explicit-directory route work independently of the 10x release.

To select a managed custom assembly, set `db.reference_db: custom_reference` and `db.custom_reference.assembly`.
When Ensembl is selected and its lookup table lacks the organism, the workflow warns before falling back to the configured custom assembly.
See [Reduced custom references](#reduced-custom-references).

For managed Cell Ranger builds, including Ensembl and custom assemblies, the workflow checksums a cached FASTA/GTF pair during planning.
With an empty cache it plans retrieval, building, counting and reporting together. The builder records the checksums once the files arrive.
The workdir pins the completed index for count. Cell Ranger reporting copies the verified original FASTA and GTF from that same index’s `gcf-source/`.
The gene table and exon intermediates live under the workdir pin’s `annotations/<source-digest>/` directory, reached through its `annotation` symlink.
A resumed run reuses the initial pin while the source pair matches. Replacing either source selects another generation even if its timestamp is preserved.

In a project using both Cell Ranger and STARsolo, STARsolo reads `<REF_DIR>/anno/genes.tsv` from its managed provider directory.
With an explicit `ref_dir`, both quantifiers read the workdir gene table derived from the checksum-verified project FASTA and GTF.
STARsolo builds its read-length-specific index under that snapshot’s `index/genome/` directory using the same FASTA and GTF.
Replacing either file selects a new index and reruns STARsolo counting and AnnData conversion, even when the replacement keeps the original timestamp.
Cell Ranger index verification runs only for targets that consume its pin, so STARsolo-only targets can run independently.

### Download dates in logs

Reference rules write their source and download date to the rule's log (e.g.
`Mirbase mature,<release>,<url>,<dd-mm-yyyy>`). The date is taken when the job
runs, in the shell command (`$(date +%d-%m-%Y)`), never in `params`: Snakemake
reruns a finished job when its params change, so a date param made completed
references and everything downstream rerun when a workflow was resumed on a later
day ([#330](https://github.com/gcfntnu/gcf-workflows/issues/330)).
`.tests/test_rule_param_dates.py` checks that no rule's `params` reads the clock.

### Langmead prebuilt Kraken databases

The shared default is full **`k2_pluspf`, release `'20260626'`**. The
[upstream catalogue](https://benlangmead.github.io/aws-indexes/k2) publishes both
that archive and the optional `k2_pluspf_16_GB` variant for the same release.
The capped variant is an explicit override, not the production default:

```yaml
db:
  langmead:
    release: '20260626'
    assembly: k2_pluspf_16_GB
```

[`langmead.db`](langmead.db) reads the effective `db.langmead.release` and
`db.langmead.assembly` values. For an assembly beginning with `k2`, it requests:

```text
https://genome-idx.s3.amazonaws.com/kraken/<assembly>_<release>.tar.gz
```

The destination is:

```text
<EXT_DIR>/langmead/release-<release>/metagenome/<assembly>/
```

Use a quoted string for the release. Select an assembly/release combination
that actually exists upstream; the getter does not invent a reduced database or
fall back to a different release. A size-limited prebuilt variant is selected
through the assembly value using the same getter.

The `k2` rule requires `hash.k2d`, `taxo.k2d`, `opts.k2d`, `seqid2taxid.map` and
Bracken distributions for read lengths 50, 75, 100, 150, 200, 250 and 300.
An archive that lacks those declared products does not satisfy the rule.
Assemblies not beginning with `k2` follow the separate legacy SILVA branch.

For the two PlusPF 20260626 variants, acquisition checks every extracted file
against the upstream per-file MD5 manifest before reporting success. The
reviewed manifests are retained in [checksums](checksums); their sources are
[full PlusPF](https://genome-idx.s3.amazonaws.com/kraken/pluspf_20260626/pluspf.md5)
and [PlusPF-16](https://genome-idx.s3.amazonaws.com/kraken/pluspf_16_GB_20260626/pluspf_16_GB.md5).
The archive entry is omitted because acquisition streams the compressed archive
directly to extraction; all 17 extracted files are checked, including the three
Kraken files, sequence mapping and all seven Bracken distributions. Existing
older/custom release acquisition remains supported without requiring these new
manifests. Existing completed cache entries are not automatically re-downloaded;
to audit them, run `md5sum --check /absolute/path/to/the/matching/manifest.md5`
from the database directory.

RAM copies now use
`/dev/shm/langmead/release-<release>/metagenome/<assembly>/`. Shared Kraken screening
and metagenome classification with a getter-managed database obtain that location from the producer.
With a project `ref_dir`, metagenome resolves the absolute database path during planning and passes it to both Kraken and Bracken.
For the required index files and read-length distribution, see the [metagenome configuration](../metagenome/README.md#key-configuration).

Old disk and RAM paths
are left intact; do not remove them while retained or active workflows need
them. New workdirs use the new default; explicit overrides and retained workflow
copies keep their configured release. RAM capacity must accommodate the selected
database (the full June 2026 index is about 110 GiB); the limited variant's size
cap changes classification sensitivity independently of the release.

Metagenome Bracken selects the nearest available distribution to the configured
read length; the existing lower-length tie behavior remains unchanged. The
available lengths are 50, 75, 100, 150, 200, 250 and 300. Focused dry-run tests
exercise both consumers, both PlusPF variants, the old release, every supported
distribution and regeneration when a required distribution is missing:

```bash
python -m unittest discover -s metagenome/.tests -p 'test_*.py' -v
```

These tests require Snakemake and pandas. They do not download references or
classify reads; scientific integration results must be recorded separately.

Check library-prep overrides as well as `main.config`: a selected kit can
override the main database selection, and explicit project configuration can
take precedence over defaults. Inspect the actual effective configuration used
by the analysis. Preserve unrelated database choices required by other routes.

A reduced database changes classification results. Test useful retained
classification and abundance outputs; do not require numerical equality with a
different database/release. Record cold download/preparation time separately
from warm workflow execution.

### Reduced custom references

The custom-reference route uses this configuration shape (names are illustrative):

```yaml
db:
  reference_db: custom_reference
  custom_reference:
    assembly: ExampleAssembly_reduced
```

Together with the normalized organism name, the FASTA location is:

```text
<EXT_DIR>/custom_reference/<organism>/<assembly>/fasta/genome.fa
```

[`custom_reference.db`](custom_reference.db) obtains missing resources from the
configured `winecellar.url` under the matching `custom_reference` path. It does
not perform downsampling. A private fixture preparation step can place a
verified derived FASTA at the selected resource path before the normal
[index rules](indexes.smk) are invoked. Required annotations remain separate
resources; a FASTA-only fixture is not sufficient for routes that need them.

For a reduced WGS fixture, retain intact reference sequences to preserve names
and coordinates, and select corresponding paired reads. Keep a deterministic
recipe with original/derived checksums and paired-read integrity checks. Build
indexes under the new assembly's paths using the existing rules. Validate
meaningful alignment and coverage rather than completion alone.
