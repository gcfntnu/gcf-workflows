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

The project reference directory can have any name and must contain the augmented FASTA and its matching annotation. Relative paths start at the project working directory.

| File | Contents |
| --- | --- |
| `refs/panel-v1/fasta/genome.fa` | Genome sequences including the reporter contig |
| `refs/panel-v1/anno/genes.gtf` | Gene annotation including the reporter feature with the same contig name and coordinates |

For a single-exon reporter spanning a 720-base FASTA contig named `reporter`, this tab-delimited GTF record supplies the matching annotation.

```text
reporter	custom	exon	1	720	.	+	.	gene_id "REPORTER"; transcript_id "REPORTER-T1";
```

When a gene has only exon records, the report derives one row in `genes.tsv` from its exons and joins it to the AnnData feature with the same `gene_id`.
Gene IDs are matched as text, preserving literal identifiers such as `NA` and `001`. Supplied gene records keep their coordinates and attributes.

The builder copies the selected FASTA and GTF into a private working directory, verifies their checksums and runs
`cellranger mkref --fasta refs/panel-v1/fasta/genome.fa --genes refs/panel-v1/anno/genes.gtf` there.
Absolute project paths and paths containing `..` are mapped into that private directory. The published index retains the verified originals in `gcf-source/`.
Count passes its workdir’s `index` symlink to `--transcriptome`. The link pins the published generation in `refs/panel-v1/cellranger/builds/<source-digest>/<build-key>`.
When `ref_dir` is a symlink, the workflow resolves the selected directory before building. Changing the symlink later cannot redirect a planned count.
A second project may select a published generation’s `gcf-source/` as `ref_dir`, directly or through a symlink.
Before planning derived jobs, the importer compares the SHA-256 checksum of each frozen source file with the owner’s `gcf-mkref.json`.
A missing or invalid manifest, or a checksum mismatch, stops the import and identifies the affected path.
In a single-cell importer, the verified FASTA and GTF are staged under `<QUANT_INTERIM>/reference/imported/<source-path-key>/<source-digest>/`.
Split-Pipe builds its index under `index/genome/splitpipe/` there and uses it for counting. The report derives its exon files and gene table from the same staged pair.
Request derived targets such as `anno/exons.gtf` under the importing workdir’s staged reference. Planning rejects writes inside a published generation.
For Cell Ranger, a matching recipe reuses and verifies the existing generation. A different genome identifier or image builds under the importing workdir’s `QUANT_INTERIM/reference/cellranger/builds/`.
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
On a resumed report, the workflow checks the workdir FASTA and GTF against the selected source checksums.
If either staged file has changed or is missing from a previously used copy, the workflow stages both files under `<staging-directory>.repairs/<n>/` and rebuilds the exon files and gene table there.
Later reports reuse that copy and the existing Cell Ranger count and index.
The configured provider and annotation settings still apply.

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
With an empty cache it plans retrieval, building, counting and reporting together. Each managed FASTA or GTF download is staged in a temporary file beside its destination and published by atomic rename once complete.
The builder records the checksums once the files arrive.
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

### miRBase GFF3 correction

[`mirbase.db`](mirbase.db) downloads the hairpin and mature sequences and the genome GFF3 of `db.mirbase.release`. The
GFF3 that every consumer reads, `EXT_DIR/mirbase/<organism>/release-<release>/anno/genes.<id>.gff3` (named after the
corrections, below), is miRBase's file
(`sources/<code>.gff3`) corrected with the reviewed rows of
[`tables/mirbase_gff_corrections.tsv`](tables/mirbase_gff_corrections.tsv) and checked by
[`mirbase_gff_correct.py`](scripts/mirbase_gff_correct.py) ([#429](https://github.com/gcfntnu/gcf-workflows/issues/429)).
Its consumers are UNITAS's annotation and `genes.<id>.tsv`, the genome routes, smrnaseq's mirtop and biocore_miraligner.

- **Table rows:** each names a GFF3 row by ID, `Derives_from` and old position, and gives its new position, how it is
  checked and the reason. A row not found at that position stops the build, so a changed download or a stale table is
  noticed. So does a row moved onto another row of the same miRNA and hairpin, which would merge two loci.
- **Check:** every row must agree with the release's sequences, or the build stops:
  - a hairpin row has a strand and spans its hairpin sequence;
  - a mature row lies where its sequence is in its hairpin;
  - chromosome names use one style.

  Two table levels relax this. `placement`, for a mature whose sequence differs from a paralogous hairpin, checks only
  that the row lies in the hairpin, on its strand, with the sequence's length. `none`, for a reviewed anomaly, skips
  the check.
- **Record:** the applied rows go to `logs/gff_corrections.<id>.tsv`.
- **Name:** `<id>` (`MIRBASE_GFF_ID`) is a digest of the table and `mirbase_gff_correct.py`, so any change to either
  gives new files, which an `EXT_DIR` then builds, and the consumers rebuild through the new path. Unchanged files are
  never rebuilt: the table is an `ancient()` input, since a checkout gives it the checkout time
  ([#438](https://github.com/gcfntnu/gcf-workflows/issues/438)). The shared files built from the GFF3 carry the same
  identity ([#440](https://github.com/gcfntnu/gcf-workflows/issues/440)): `anno/genes.<id>.tsv`, and a `gff-<id>`
  folder in the mirflowz (chromosome map, extended annotations), tiny-count, TEsmall and biocore_miraligner
  references. Code versions with different corrections can therefore share an `EXT_DIR`, each using its own; the
  files built from the genome alone (the genome FASTA and its indexes, mirflowz's segemehl indexes) are shared.

**miRBase 23.1's corrections:**
- hairpin rows that cover only the stem, with their matures placed from them; some of these rows have no strand;
- hairpin rows not updated when release 23 changed the hairpin sequence;
- mature rows away from their sequence;
- in mouse, rows still on GRCm38 instead of GRCm39.

The new positions are where the hairpin sequences lie exactly in GRCh38, GRCm39 and Rnor_6.0. Matures whose sequence
differs from a paralogous hairpin are placed at their best match and checked for placement only. Releases 20, 21 and 22
pass unchanged; 22 checks three rat rows for placement only, as their accession is missing from its `mature.fa`.

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
  of it, and once per UNITAS dump, from which the piRNA clusters come. The genome FASTA and its index are built from the
  genome alone, so a new annotation does not rebuild them or the biocore_htseq index
  ([#433](https://github.com/gcfntnu/gcf-workflows/issues/433)). With `quant.tesmall.structural_rnacentral`, a database
  derived from that one adds RNAcentral's snoRNA, scaRNA, Y RNA, vault RNA and 7SK entries (release 28.0) to the
  structural RNA class, in its own folder beside the base database, whose other files it links
  ([`tesmall_rnacentral.py`](scripts/tesmall_rnacentral.py), [#378](https://github.com/gcfntnu/gcf-workflows/issues/378));
  human and mouse only.
- **tiny-count reference** ([`tinyrna.db`](tinyrna.db), [`tinyrna_gff.py`](scripts/tinyrna_gff.py)): TEsmall's Bowtie 1
  genome index, and the miRBase GFF renamed to the genome's sequence names with the mirflowz chromosome map
  (`mirna.v<n>.gff3`; `<n>` is raised when the script's output changes, so an older GFF is rebuilt). A mature miRNA
  listed at one locus under several hairpins is one row with every parent; repeated IDs at other loci get a numeric
  suffix.
- **biocore_htseq reference** ([`biocore_htseq.db`](biocore_htseq.db)): a Bowtie 2 index of TEsmall's genome FASTA, under
  `EXT_DIR/biocore_htseq/<mirflowz reference>/genome/bowtie2/`; the route annotates with tiny-count's renamed miRBase GFF.
- **biocore_miraligner database** ([`biocore_miraligner.db`](biocore_miraligner.db), [`miraligner_db.py`](scripts/miraligner_db.py)): the
  organism's miRBase hairpins and a `miRNA.str` with each mature miRNA's position in its hairpin, where its sequence lies
  there, guided by the GFF3 (for miRBase 22 identical to miRBase's own file; miRBase serves none for 23.1). The record
  `miraligner_db.tsv` lists the rows whose GFF3 position disagrees with the sequences; with the corrected GFF3 (below),
  those the correction checks for placement only or not at all.

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
