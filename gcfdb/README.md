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

### Download dates in logs

Reference rules write their source and download date to the rule's log (e.g.
`Mirbase mature,<release>,<url>,<dd-mm-yyyy>`). The date is taken when the job
runs, in the shell command (`$(date +%d-%m-%Y)`), never in `params`: Snakemake
reruns a finished job when its params change, so a date param made completed
references and everything downstream rerun when a workflow was resumed on a later
day ([#330](https://github.com/gcfntnu/gcf-workflows/issues/330)).
`.tests/test_rule_param_dates.py` checks that no rule's `params` reads the clock.

### Small-RNA quantifier references

The smallrna quantification routes each use their own reference under `EXT_DIR`; the getter files' docstrings give
their identity, sources and checks.

- **UNITAS reference dumps** ([`unitas.db`](unitas.db)): `UNITAS_refdump_<organism>` from `db.unitas.release`
  (miRBase) or `db.unitas.mirgenedb_release` (MirGeneDB, with `db.unitas.mirna_db: mirgenedb`), for human, mouse and
  rat only. `unitas_db_verify` checks a cached dump against its release's `SHA256SUMS`. smallrna's `unitas_reference`
  (per project, [`unitas_reference.py`](../smallrna/rules/quant/scripts/unitas_reference.py)) records the dump's
  releases and stops when its miRNA database or miRBase release differs from `db.unitas.mirna_db` or
  `db.mirbase.release`: UNITAS counts are joined to the annotation by name, and names change between releases (with
  the miRBase 22 GFF, 35% of the miRNAs a miRBase 23.1 dump counts are unannotated).
- **miRge3 libraries** ([`mirge3.db`](mirge3.db), [`mirge3_library.py`](scripts/mirge3_library.py)): dated releases
  built from public sources ([gcfntnu/gcf-docker#41](https://github.com/gcfntnu/gcf-docker/pull/41), which lists the
  changes from the old SourceForge libraries: miRBase 23.1, MirGeneDB 3.0, Ensembl 116, merge groups of identical
  sequences only), checked against `SHA256SUMS`; `sourceforge` keeps the old unversioned library for retained
  analyses. smallrna builds each project's own library with its spike-in index
  ([`mirge3_project_library.py`](../smallrna/rules/quant/scripts/mirge3_project_library.py)).
- **mirflowz references** ([`mirflowz.db`](mirflowz.db)): genome, GTF and miRBase GFF on one assembly, the Ensembl
  release fixed per organism in `db.mirflowz.ensembl_release`, prepared with mirflowz's own rules. segemehl needs about
  60 GB of memory for a human or mouse genome index.
- **TEsmall databases** ([`tesmall.db`](tesmall.db)): built on the mirflowz reference's genome, with a Bowtie 1 index
  of it, and once per UNITAS dump, from which the piRNA clusters come.
- **tiny-count reference** ([`tinyrna.db`](tinyrna.db), [`tinyrna_gff.py`](scripts/tinyrna_gff.py)): TEsmall's Bowtie 1
  genome index, and the miRBase GFF renamed to the genome's sequence names with the mirflowz chromosome map
  (`mirna.v<n>.gff3`; `<n>` is raised when the script's output changes, so an older GFF is rebuilt). A mature miRNA
  listed at one locus under several hairpins is one row with every parent; repeated IDs at other loci get a numeric
  suffix.
- **miraligner database** ([`miraligner.db`](miraligner.db), [`miraligner_db.py`](scripts/miraligner_db.py)): the
  organism's miRBase hairpins and a `miRNA.str` with each mature miRNA's position in its hairpin, where its sequence lies
  there, guided by the GFF3 (for miRBase 22 identical to miRBase's own file; miRBase serves none for 23.1). The record
  `miraligner_db.tsv` lists the miRBase 23.1 rows whose GFF3 position disagrees with the sequences.

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
`/dev/shm/langmead/release-<release>/metagenome/<assembly>/`. Both shared Kraken
screening and metagenome classification obtain that location from the producer,
so neither can silently use another release's RAM copy. Old disk and RAM paths
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
