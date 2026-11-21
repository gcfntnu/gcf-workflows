#-*- mode:snakemake -*-
from functools import lru_cache
import stat

include: 'scripts/reference_identity.py'


@lru_cache(maxsize=4096)
def has_publication_marker(directory, directory_version):
    # Cache marker checks while the directory entries are unchanged, including between
    # DAG construction and execution after a dependent job publishes a new generation.
    if ((directory / 'gcf-source').is_dir() and (directory / 'reference.json').is_file()
            and (directory / 'star/SA').is_file()):
        return True
    # A copied manifest alone may just be project provenance. A relocated generation
    # also has at least one of its published source or index members.
    if ((directory / CELLRANGER_MANIFEST).is_file()
            and any((directory / member).is_file() for member in ('gcf-source/fasta/genome.fa', 'gcf-source/anno/genes.gtf', 'star/SA'))):
        try:
            read_cellranger_manifest(directory)
        except (FileNotFoundError, ValueError):
            return False
        return True
    return False


def is_publication_directory(path):
    if is_managed_cellranger_generation(path):
        return True
    try:
        info = path.stat()
    except (FileNotFoundError, NotADirectoryError):
        return False
    return stat.S_ISDIR(info.st_mode) and has_publication_marker(path, (info.st_dev, info.st_ino, info.st_ctime_ns))


# Each job owns its output directory entries, even when another project published the reference.
workdir_reference_links = set()
workdir_reference_roots = {}


def check_reference_outputs(wildcards, rulename):
    # Canonicalize each output directory once per job. Re-evaluate symlinks at execution
    # because dependencies can create a new pin after the initial plan.
    @lru_cache(maxsize=None)
    def published_ancestor(path):
        if is_publication_directory(path):
            return path
        return published_ancestor(path.parent) if path != path.parent else None

    resolve_parent = lru_cache(maxsize=None)(Path.resolve)
    outputs, _ = workflow.get_rule(rulename).expand_output(wildcards)
    for output in outputs:
        path = Path(output)
        entry = resolve_parent(path.parent) / path.name
        target = entry
        generation = published_ancestor(entry.parent)
        if not generation and (rulename, entry) not in workdir_reference_links:
            target = entry.resolve() if entry.is_symlink() else entry
            generation = published_ancestor(target)
        if generation:
            if generation in workdir_reference_roots:
                relative = target.relative_to(generation)
                if relative.parts[:1] == ('gcf-source',):
                    relative = relative.relative_to('gcf-source')
                advice = f'Request {workdir_reference_roots[generation] / relative} instead.'
            else:
                advice = f'Select {generation / "gcf-source"} as ref_dir to build derived files in the importing workdir.'
            raise ValueError(f'Cannot write inside a published Cell Ranger generation: {output}. {advice}')
    return 0


# Snakemake evaluates default resources before preparing outputs in execution and dry runs.
# Keep any existing rule, profile or command-line resource with the same name intact.
guard_resource = '_gcf_cellranger_output_guard'
resource_settings = workflow.resource_settings
used_resources = set(resource_settings.default_resources.parsed) | set(resource_settings.resources)
for configured in resource_settings.overwrite_resources.values():
    used_resources.update(configured)
for registered_rule in workflow.rules:
    used_resources.update(registered_rule.resources)
while guard_resource in used_resources:
    guard_resource += '_'
resource_settings.default_resources.parsed[guard_resource] = check_reference_outputs
for registered_rule in workflow.rules:
    registered_rule.resources[guard_resource] = check_reference_outputs
