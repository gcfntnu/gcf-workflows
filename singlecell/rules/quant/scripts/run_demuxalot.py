#!/usr/bin/env python
"""Run Demuxalot genotype demultiplexing.

Supported genotype modes
------------------------
fixed
    Use the supplied donor genotypes without refinement and call
    Demultiplexer.predict_posteriors().

learn
    Count SNPs at positions present in the supplied donor VCF and refine
    donor genotypes with Demultiplexer.learn_genotypes().

learn-new-snps
    Detect additional informative SNP positions from the single-cell BAM,
    add them as weak genotype priors, recount SNPs using the expanded
    position set, and refine genotypes with
    Demultiplexer.learn_genotypes().

The output follows the GCF demultiplexing interface:

    Barcode
    doublet_type
    donor_id
    best_singlet
    doublet_score

Additional Demuxalot-specific posterior diagnostics are retained.

References
----------
Demuxalot:
https://github.com/arogozhnikov/demuxalot

Examples:
1-plain_demultiplexing.py
2-with-detection-of-new-SNPs.ipynb
3-plain_demultiplexing-with-custom-tags.py
"""

import argparse
import os
from collections import Counter
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from demuxalot import (
    BarcodeHandler,
    Demultiplexer,
    ProbabilisticGenotypes,
    count_snps,
    detect_snps_positions,
)
from demuxalot.cellranger_specific import parse_read as cellranger_parse_read


GENOTYPE_MODES = (
    "fixed",
    "learn",
    "refine",
    "learn-new-snps",
)


def parse_args():
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed and validated command-line arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "-b",
        "--barcode-filename",
        required=True,
        help="Single-cell barcode file.",
    )
    parser.add_argument(
        "-i",
        "--bam-filename",
        required=True,
        help="Single-cell BAM file.",
    )
    parser.add_argument(
        "-v",
        "--vcf-filename",
        required=True,
        help="Reference donor VCF.",
    )
    parser.add_argument(
        "-g",
        "--donor-ids",
        default=None,
        help=(
            "Comma-separated donor IDs or a file containing one donor ID "
            "per line. If omitted, all samples in the VCF are used."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output droplet_type.tsv file.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=-1,
        help="Number of joblib workers. -1 uses all available CPUs.",
    )

    parser.add_argument(
        "--genotype-mode",
        choices=GENOTYPE_MODES,
        default="learn",
        help=(
            "Genotype handling mode. 'fixed' uses the supplied genotypes "
            "without refinement; 'learn' refines genotypes at known VCF "
            "positions; 'learn-new-snps' additionally discovers informative "
            "SNP positions from the single-cell BAM. Default: learn."
        ),
    )
    parser.add_argument(
        "--doublet-prior",
        type=float,
        default=0.15,
        help=(
            "Prior probability of a doublet passed to Demuxalot. "
            "Setting this to zero disables doublet hypotheses. "
            "Default: 0.15."
        ),
    )
    parser.add_argument(
        "--min-prob",
        type=float,
        default=0.9,
        help=(
            "Minimum posterior probability of the single best Demuxalot "
            "hypothesis required for a categorical assignment. "
            "Default: 0.9."
        ),
    )
    parser.add_argument(
        "--vcf-prior-strength",
        type=float,
        default=100.0,
        help=(
            "Prior strength assigned to genotypes loaded from the VCF. "
            "Default: 100."
        ),
    )

    #
    # Custom BAM tags.
    #
    parser.add_argument(
        "--cell-tag",
        default=None,
        help=(
            "Custom BAM cell-barcode tag. If omitted, Demuxalot's "
            "Cell Ranger default is used."
        ),
    )
    parser.add_argument(
        "--umi-tag",
        default=None,
        help=(
            "Custom BAM UMI tag. If omitted, Demuxalot's Cell Ranger "
            "default is used."
        ),
    )

    #
    # New-SNP discovery.
    #
    parser.add_argument(
        "--new-snp-min-coverage",
        type=int,
        default=200,
        help=(
            "Minimum coverage used when detecting new SNP positions. "
            "Default: 200."
        ),
    )
    parser.add_argument(
        "--new-snp-min-alt-coverage",
        type=int,
        default=10,
        help=(
            "Minimum alternative-allele coverage used when detecting new "
            "SNP positions. Default: 10."
        ),
    )
    parser.add_argument(
        "--new-snp-min-alt-fraction",
        type=float,
        default=0.01,
        help=(
            "Minimum alternative-allele fraction used during new-SNP "
            "detection. Default: 0.01."
        ),
    )
    parser.add_argument(
        "--new-snp-best-per-donor",
        type=int,
        default=100,
        help=(
            "Number of best newly detected SNPs retained per donor. "
            "Default: 100."
        ),
    )
    parser.add_argument(
        "--new-snp-additional-best",
        type=int,
        default=1000,
        help=(
            "Number of additional best SNPs retained across donors. "
            "Default: 1000."
        ),
    )
    parser.add_argument(
        "--new-snp-prior-strength",
        type=float,
        default=1.0,
        help=(
            "Prior strength assigned to newly detected SNP positions. "
            "Default: 1."
        ),
    )
    parser.add_argument(
        "--new-snps-output",
        default=None,
        help=(
            "Output beta-prior file for newly detected SNP positions. "
            "If omitted in learn-new-snps mode, new_snps.betas is written "
            "beside the main output."
        ),
    )

    #
    # Optional diagnostic / reusable outputs.
    #
    parser.add_argument(
        "--learned-genotypes-output",
        default=None,
        help=(
            "Optional output path for learned Demuxalot genotype betas. "
            "Only valid for learn and learn-new-snps modes."
        ),
    )
    parser.add_argument(
        "--posterior-output",
        default=None,
        help=(
            "Optional output path for the complete Demuxalot posterior "
            "probability matrix. '.gz' is supported by pandas."
        ),
    )

    args = parser.parse_args()

    if not 0.0 <= args.doublet_prior < 1.0:
        parser.error(
            "--expected-doublet-rate must satisfy 0 <= value < 1."
        )

    if not 0.0 <= args.min_prob <= 1.0:
        parser.error("--min-prob must be between 0 and 1.")

    if args.vcf_prior_strength <= 0:
        parser.error("--vcf-prior-strength must be > 0.")

    if args.new_snp_min_coverage < 1:
        parser.error("--new-snp-min-coverage must be >= 1.")

    if args.new_snp_min_alt_coverage < 1:
        parser.error("--new-snp-min-alt-coverage must be >= 1.")

    if not 0.0 <= args.new_snp_min_alt_fraction <= 1.0:
        parser.error(
            "--new-snp-min-alt-fraction must be between 0 and 1."
        )

    if args.new_snp_best_per_donor < 1:
        parser.error("--new-snp-best-per-donor must be >= 1.")

    if args.new_snp_additional_best < 0:
        parser.error("--new-snp-additional-best must be >= 0.")

    if args.new_snp_prior_strength <= 0:
        parser.error("--new-snp-prior-strength must be > 0.")

    if (
        args.genotype_mode == "fixed"
        and args.learned_genotypes_output is not None
    ):
        parser.error(
            "--learned-genotypes-output cannot be used with "
            "--genotype-mode fixed."
        )

    if (
        args.genotype_mode != "learn-new-snps"
        and args.new_snps_output is not None
    ):
        parser.error(
            "--new-snps-output requires "
            "--genotype-mode learn-new-snps."
        )

    return args


def ensure_parent_directory(filename):
    """Create the parent directory of a filename if necessary.

    Parameters
    ----------
    filename : str or pathlib.Path
        Output filename.
    """
    parent = Path(filename).parent
    parent.mkdir(parents=True, exist_ok=True)


def ensure_bam_index(bam_filename):
    """Ensure that a BAM index exists.

    Parameters
    ----------
    bam_filename : str
        BAM filename.
    """
    with pysam.AlignmentFile(bam_filename, "rb") as bam:
        has_index = bam.has_index()

    if has_index:
        return

    print(f"Indexing BAM file: {bam_filename}")
    pysam.index(bam_filename)


def get_vcf_donors(vcf_filename):
    """Return donor sample names present in a VCF.

    Parameters
    ----------
    vcf_filename : str
        Input VCF.

    Returns
    -------
    list of str
        VCF sample names.
    """
    with pysam.VariantFile(vcf_filename) as vcf:
        donors = [str(x) for x in vcf.header.samples]

    if not donors:
        raise ValueError(
            f"No genotype samples found in VCF: {vcf_filename}"
        )

    return donors


def resolve_donor_ids(vcf_filename, donor_ids_arg):
    """Resolve donor IDs and validate them against the VCF.

    Parameters
    ----------
    vcf_filename : str
        Input donor VCF.
    donor_ids_arg : str or None
        Comma-separated donor IDs, filename containing donor IDs,
        one donor ID, or None.

    Returns
    -------
    list of str
        Sorted donor IDs.
    """
    vcf_donors = get_vcf_donors(vcf_filename)

    if donor_ids_arg is None:
        donor_ids = list(vcf_donors)

    elif os.path.isfile(donor_ids_arg):
        with open(donor_ids_arg) as fh:
            donor_ids = [
                line.strip()
                for line in fh
                if line.strip()
            ]

    elif "," in donor_ids_arg:
        donor_ids = [
            donor.strip()
            for donor in donor_ids_arg.split(",")
            if donor.strip()
        ]

    else:
        donor_ids = [donor_ids_arg.strip()]

    if not donor_ids:
        raise ValueError("No donor IDs were supplied.")

    if len(donor_ids) != len(set(donor_ids)):
        raise ValueError(
            f"Duplicate donor IDs supplied: {donor_ids}"
        )

    if any("+" in donor for donor in donor_ids):
        raise ValueError(
            "Donor IDs containing '+' are not supported because "
            "Demuxalot uses '+' to represent doublet hypotheses."
        )

    missing = sorted(set(donor_ids) - set(vcf_donors))

    if missing:
        raise ValueError(
            "Requested donors are absent from the VCF: "
            + ", ".join(missing)
        )

    return sorted(donor_ids)


def make_barcode_handler(barcode_filename, cell_tag=None):
    """Create a Demuxalot barcode handler.

    Parameters
    ----------
    barcode_filename : str
        Barcode file.
    cell_tag : str or None
        Custom BAM cell-barcode tag.

    Returns
    -------
    BarcodeHandler
        Configured barcode handler.
    """
    if cell_tag is None:
        return BarcodeHandler.from_file(barcode_filename)

    return BarcodeHandler.from_file(
        barcode_filename,
        tag=cell_tag,
    )


def make_parse_read(umi_tag=None):
    """Create the BAM-read parser.

    Parameters
    ----------
    umi_tag : str or None
        Custom UMI tag.

    Returns
    -------
    callable
        Demuxalot-compatible read parser.
    """
    if umi_tag is None:
        return cellranger_parse_read

    return partial(
        cellranger_parse_read,
        umi_tag=umi_tag,
    )


def make_genotypes(vcf_filename, donor_ids, prior_strength):
    """Load donor genotypes from VCF.

    Parameters
    ----------
    vcf_filename : str
        Donor VCF.
    donor_ids : list of str
        Donors to load.
    prior_strength : float
        VCF genotype prior strength.

    Returns
    -------
    ProbabilisticGenotypes
        Loaded donor genotypes.
    """
    genotypes = ProbabilisticGenotypes(
        genotype_names=donor_ids
    )

    genotypes.add_vcf(
        vcf_filename,
        prior_strength=prior_strength,
    )

    return genotypes


def count_genotype_snps(
    bam_filename,
    genotypes,
    barcode_handler,
    parse_read,
    threads,
):
    """Count single-cell SNP evidence at genotype positions.

    Parameters
    ----------
    bam_filename : str
        Single-cell BAM.
    genotypes : ProbabilisticGenotypes
        Genotype object defining SNP positions.
    barcode_handler : BarcodeHandler
        Barcode handler.
    parse_read : callable
        BAM-read parser.
    threads : int
        Number of joblib workers.

    Returns
    -------
    dict
        Demuxalot compressed SNP counts.
    """
    return count_snps(
        bamfile_location=bam_filename,
        chromosome2positions=genotypes.get_chromosome2positions(),
        barcode_handler=barcode_handler,
        parse_read=parse_read,
        joblib_n_jobs=threads,
    )


def run_fixed(
    args,
    genotypes,
    barcode_handler,
    parse_read,
):
    """Run Demuxalot without genotype refinement.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments.
    genotypes : ProbabilisticGenotypes
        Input genotypes.
    barcode_handler : BarcodeHandler
        Barcode handler.
    parse_read : callable
        Read parser.

    Returns
    -------
    ProbabilisticGenotypes
        Input genotypes.
    pandas.DataFrame
        Posterior probabilities.
    dict
        SNP counts.
    """
    snps = count_genotype_snps(
        bam_filename=args.bam_filename,
        genotypes=genotypes,
        barcode_handler=barcode_handler,
        parse_read=parse_read,
        threads=args.threads,
    )

    _, posterior_probabilities = (
        Demultiplexer.predict_posteriors(
            snps,
            genotypes=genotypes,
            barcode_handler=barcode_handler,
            doublet_prior=args.doublet_prior,
        )
    )

    return genotypes, posterior_probabilities, snps


def run_learn(
    args,
    genotypes,
    barcode_handler,
    parse_read,
):
    """Run Demuxalot with genotype refinement at known VCF SNPs.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments.
    genotypes : ProbabilisticGenotypes
        Input genotypes.
    barcode_handler : BarcodeHandler
        Barcode handler.
    parse_read : callable
        Read parser.

    Returns
    -------
    ProbabilisticGenotypes
        Learned genotypes.
    pandas.DataFrame
        Posterior probabilities returned by learn_genotypes().
    dict
        SNP counts.
    """
    snps = count_genotype_snps(
        bam_filename=args.bam_filename,
        genotypes=genotypes,
        barcode_handler=barcode_handler,
        parse_read=parse_read,
        threads=args.threads,
    )

    learned_genotypes, posterior_probabilities = (
        Demultiplexer.learn_genotypes(
            snps,
            genotypes=genotypes,
            barcode_handler=barcode_handler,
            doublet_prior=args.doublet_prior,
        )
    )

    return learned_genotypes, posterior_probabilities, snps

def run_refine(
    args,
    genotypes,
    barcode_handler,
    parse_read,
):
    """Refine donor genotypes as singlets before modelling doublets.

    Genotypes are first learned with doublet hypotheses disabled. The
    refined genotypes are then used for final posterior prediction with
    the configured doublet prior.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments.
    genotypes : ProbabilisticGenotypes
        Initial donor genotypes.
    barcode_handler : BarcodeHandler
        Barcode handler.
    parse_read : callable
        Read parser.

    Returns
    -------
    ProbabilisticGenotypes
        Refined donor genotypes.
    pandas.DataFrame
        Final posterior probabilities.
    dict
        SNP counts.
    """
    snps = count_genotype_snps(
        bam_filename=args.bam_filename,
        genotypes=genotypes,
        barcode_handler=barcode_handler,
        parse_read=parse_read,
        threads=args.threads,
    )

    print(
        "Refining donor genotypes with doublet hypotheses disabled ..."
    )

    refined_genotypes, _ = Demultiplexer.learn_genotypes(
        snps,
        genotypes=genotypes,
        barcode_handler=barcode_handler,
        doublet_prior=0.0,
    )

    print(
        "Predicting final posteriors with "
        f"doublet_prior={args.doublet_prior} ..."
    )

    _, posterior_probabilities = Demultiplexer.predict_posteriors(
        snps,
        genotypes=refined_genotypes,
        barcode_handler=barcode_handler,
        doublet_prior=args.doublet_prior,
    )

    return refined_genotypes, posterior_probabilities, snps

def run_learn_new_snps(
    args,
    genotypes,
    barcode_handler,
    parse_read,
    output_dir,
):
    """Run genotype refinement with additional SNP discovery.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments.
    genotypes : ProbabilisticGenotypes
        Input genotypes.
    barcode_handler : BarcodeHandler
        Barcode handler.
    parse_read : callable
        Read parser.
    output_dir : pathlib.Path
        Main output directory.

    Returns
    -------
    ProbabilisticGenotypes
        Learned genotypes.
    pandas.DataFrame
        Posterior probabilities returned by learn_genotypes().
    dict
        SNP counts using known plus newly detected SNP positions.
    """
    if args.new_snps_output is None:
        new_snps_filename = output_dir / "new_snps.betas"
    else:
        new_snps_filename = Path(args.new_snps_output)

    ensure_parent_directory(new_snps_filename)

    print("Detecting additional informative SNP positions ...")

    selected_snps = detect_snps_positions(
        bamfile_location=args.bam_filename,
        genotypes=genotypes,
        barcode_handler=barcode_handler,
        minimum_coverage=args.new_snp_min_coverage,
        minimum_alternative_fraction=(
            args.new_snp_min_alt_fraction
        ),
        minimum_alternative_coverage=(
            args.new_snp_min_alt_coverage
        ),
        n_best_snps_per_donor=(
            args.new_snp_best_per_donor
        ),
        n_additional_best_snps=(
            args.new_snp_additional_best
        ),
        parse_read=parse_read,
        joblib_n_jobs=args.threads,
        result_beta_prior_filename=str(
            new_snps_filename
        ),
    )

    if not selected_snps:
        raise RuntimeError(
            "Demuxalot did not detect any additional SNP positions. "
            "Use --genotype-mode learn if new-SNP discovery is not "
            "appropriate for this dataset."
        )

    print(
        f"Detected {len(selected_snps)} additional SNP positions."
    )
    print(
        f"New-SNP beta priors: {new_snps_filename}"
    )

    enriched_genotypes = genotypes.clone()

    enriched_genotypes.add_prior_betas(
        str(new_snps_filename),
        prior_strength=args.new_snp_prior_strength,
    )

    snps = count_genotype_snps(
        bam_filename=args.bam_filename,
        genotypes=enriched_genotypes,
        barcode_handler=barcode_handler,
        parse_read=parse_read,
        threads=args.threads,
    )

    learned_genotypes, posterior_probabilities = (
        Demultiplexer.learn_genotypes(
            snps,
            genotypes=enriched_genotypes,
            barcode_handler=barcode_handler,
            doublet_prior=args.doublet_prior,
        )
    )

    return learned_genotypes, posterior_probabilities, snps


def validate_posterior_probabilities(
    posterior_probabilities,
    donor_ids,
):
    """Validate Demuxalot posterior probabilities.

    The posterior matrix is expected to contain one column for every
    requested singlet donor and, when doublet modelling is enabled,
    additional donor-pair hypotheses. Posterior probabilities must sum
    to approximately one for every barcode.

    Parameters
    ----------
    posterior_probabilities : pandas.DataFrame
        Demuxalot posterior probability matrix.
    donor_ids : list of str
        Expected singlet donor names.

    Returns
    -------
    pandas.DataFrame
        Validated posterior matrix.
    list of str
        Singlet columns.
    list of str
        Doublet columns.
    """
    if not isinstance(
        posterior_probabilities,
        pd.DataFrame,
    ):
        raise TypeError(
            "Demuxalot posterior probabilities must be a "
            "pandas DataFrame."
        )

    if posterior_probabilities.empty:
        raise ValueError(
            "Demuxalot returned an empty posterior matrix."
        )

    posterior_probabilities = (
        posterior_probabilities.copy()
    )

    posterior_probabilities.columns = (
        posterior_probabilities.columns.map(str)
    )
    posterior_probabilities.index = (
        posterior_probabilities.index.map(str)
    )

    if posterior_probabilities.columns.duplicated().any():
        duplicates = (
            posterior_probabilities.columns[
                posterior_probabilities.columns.duplicated()
            ]
            .unique()
            .tolist()
        )
        raise ValueError(
            "Duplicate posterior columns returned by Demuxalot: "
            f"{duplicates}"
        )

    missing_donors = [
        donor
        for donor in donor_ids
        if donor not in posterior_probabilities.columns
    ]

    if missing_donors:
        raise ValueError(
            "Expected donor columns are absent from the Demuxalot "
            "posterior matrix: "
            + ", ".join(missing_donors)
        )

    singlet_columns = list(donor_ids)

    doublet_columns = [
        column
        for column in posterior_probabilities.columns
        if column not in singlet_columns
    ]

    malformed_doublet_columns = [
        column
        for column in doublet_columns
        if "+" not in column
    ]

    if malformed_doublet_columns:
        raise ValueError(
            "Unexpected non-singlet posterior columns returned by "
            "Demuxalot: "
            + ", ".join(malformed_doublet_columns)
        )

    values = posterior_probabilities.to_numpy(
        dtype=float
    )

    if not np.isfinite(values).all():
        raise ValueError(
            "Demuxalot posterior matrix contains non-finite values."
        )

    tolerance = 1e-8

    if values.min() < -tolerance:
        raise ValueError(
            "Demuxalot posterior matrix contains negative "
            "probabilities."
        )

    if values.max() > 1.0 + tolerance:
        raise ValueError(
            "Demuxalot posterior matrix contains probabilities > 1."
        )

    posterior_sum = posterior_probabilities.sum(axis=1)

    if not np.allclose(
        posterior_sum.to_numpy(),
        1.0,
        rtol=1e-6,
        atol=1e-6,
    ):
        deviation = (
            posterior_sum - 1.0
        ).abs()

        raise ValueError(
            "Demuxalot posterior probabilities do not sum to 1 "
            "per barcode. Maximum absolute deviation: "
            f"{deviation.max():.6g}"
        )

    return (
        posterior_probabilities,
        singlet_columns,
        doublet_columns,
    )


def make_droplet_table(
    posterior_probabilities,
    donor_ids,
    min_prob,
):
    """Convert Demuxalot posteriors to the GCF demultiplexing format.

    Categorical assignments follow the Demuxalot example logic:
    the single best hypothesis must exceed ``min_prob``.

    ``doublet_score`` is the total posterior probability assigned to
    donor-pair hypotheses. It is therefore continuous and independent
    of whether one particular donor pair passes the categorical
    assignment threshold.

    Parameters
    ----------
    posterior_probabilities : pandas.DataFrame
        Demuxalot posterior probability matrix.
    donor_ids : list of str
        Singlet donor names.
    min_prob : float
        Minimum posterior probability for categorical assignment.

    Returns
    -------
    pandas.DataFrame
        Per-barcode demultiplexing results.
    """
    (
        posterior_probabilities,
        singlet_columns,
        doublet_columns,
    ) = validate_posterior_probabilities(
        posterior_probabilities,
        donor_ids,
    )

    singlet_probabilities = posterior_probabilities[
        singlet_columns
    ]

    best_singlet = singlet_probabilities.idxmax(
        axis=1
    )
    best_singlet_prob = singlet_probabilities.max(
        axis=1
    )

    best_assignment = posterior_probabilities.idxmax(
        axis=1
    )
    best_assignment_prob = posterior_probabilities.max(
        axis=1
    )

    if doublet_columns:
        doublet_probabilities = posterior_probabilities[
            doublet_columns
        ]

        doublet_score = doublet_probabilities.sum(
            axis=1
        )
        best_doublet = doublet_probabilities.idxmax(
            axis=1
        )
        best_doublet_prob = doublet_probabilities.max(
            axis=1
        )

    else:
        doublet_score = pd.Series(
            0.0,
            index=posterior_probabilities.index,
            dtype=float,
        )
        best_doublet = pd.Series(
            pd.NA,
            index=posterior_probabilities.index,
            dtype="object",
        )
        best_doublet_prob = pd.Series(
            0.0,
            index=posterior_probabilities.index,
            dtype=float,
        )

    #
    # The upstream Demuxalot example uses:
    #
    # probs[probs.max(axis=1).gt(0.9)].idxmax(axis=1)
    #
    # Keep the strict ">" semantics here.
    #
    confident = best_assignment_prob > min_prob

    best_is_singlet = best_assignment.isin(
        singlet_columns
    )

    singlet_call = confident & best_is_singlet
    doublet_call = confident & ~best_is_singlet

    result = pd.DataFrame(
        index=posterior_probabilities.index
    )
    result.index.name = "Barcode"

    result["doublet_type"] = "unassigned"
    result["donor_id"] = "unassigned"

    result.loc[
        singlet_call,
        "doublet_type",
    ] = "singlet"

    result.loc[
        doublet_call,
        "doublet_type",
    ] = "doublet"

    result.loc[
        singlet_call,
        "donor_id",
    ] = best_assignment.loc[
        singlet_call
    ].astype(str)

    result.loc[
        doublet_call,
        "donor_id",
    ] = "doublet"

    #
    # Common GCF demultiplexing interface.
    #
    result["best_singlet"] = best_singlet.astype(
        str
    )
    result["doublet_score"] = doublet_score.astype(
        float
    )

    #
    # Demuxalot-specific diagnostic information.
    #
    result["best_singlet_prob"] = (
        best_singlet_prob.astype(float)
    )
    result["best_assignment"] = (
        best_assignment.astype(str)
    )
    result["best_assignment_prob"] = (
        best_assignment_prob.astype(float)
    )
    result["best_doublet"] = best_doublet
    result["best_doublet_prob"] = (
        best_doublet_prob.astype(float)
    )

    return result


def summarize_counted_snps(snp_counts):
    """Summarize SNP evidence and generate coverage histograms.

    Parameters
    ----------
    snp_counts : dict
        SNP counts returned by Demuxalot.

    Returns
    -------
    pandas.DataFrame
        Per-chromosome SNP summary.
    matplotlib.figure.Figure
        SNP-call and transcript-count histograms.
    """
    from matplotlib import pyplot as plt

    records = []

    barcode2number_of_calls = Counter()
    barcode2number_of_transcripts = Counter()

    for chromosome, calls in snp_counts.items():
        records.append(
            {
                "chromosome": chromosome,
                "n_molecules": calls.n_molecules,
                "n_snp_calls": calls.n_snp_calls,
            }
        )

        barcode2number_of_transcripts.update(
            Counter(
                calls.molecules["compressed_cb"]
            )
        )

        barcodes = calls.molecules[
            "compressed_cb"
        ][
            calls.snp_calls["molecule_index"]
        ]

        barcode2number_of_calls.update(
            Counter(barcodes)
        )

    if not records:
        raise ValueError(
            "Demuxalot did not return SNP counts for any chromosome."
        )

    def truncate_at_percentile(
        values,
        percentile=99.5,
    ):
        """Clip values for histogram visualization."""
        values = np.asarray(
            list(values),
            dtype=float,
        )

        if values.size == 0:
            return values

        return values.clip(
            0,
            np.percentile(
                values,
                percentile,
            ),
        )

    fig, (ax1, ax2) = plt.subplots(
        ncols=2,
        figsize=(12, 5),
    )

    ax1.hist(
        truncate_at_percentile(
            barcode2number_of_calls.values()
        ),
        histtype="step",
        bins=20,
    )
    ax1.set_ylabel("barcodes")
    ax1.set_xlabel("SNP calls per droplet")

    ax2.hist(
        truncate_at_percentile(
            barcode2number_of_transcripts.values()
        ),
        histtype="step",
        bins=20,
    )
    ax2.set_ylabel("number of barcodes")
    ax2.set_xlabel("transcripts per droplet")

    fig.tight_layout()

    summary = (
        pd.DataFrame(records)
        .sort_values("chromosome")
        .set_index("chromosome")
    )

    return summary, fig


def run_demux(args):
    """Run the selected Demuxalot workflow.

    Parameters
    ----------
    args : argparse.Namespace
        Validated command-line arguments.

    Returns
    -------
    pandas.DataFrame
        Per-barcode demultiplexing results.
    pandas.DataFrame
        Per-chromosome SNP summary.
    matplotlib.figure.Figure
        SNP coverage figure.
    ProbabilisticGenotypes
        Final genotype model.
    pandas.DataFrame
        Full posterior probability matrix.
    """
    donor_ids = resolve_donor_ids(
        args.vcf_filename,
        args.donor_ids,
    )

    print(
        "Demuxalot donors: "
        + ", ".join(donor_ids)
    )
    print(
        f"Genotype mode: {args.genotype_mode}"
    )

    barcode_handler = make_barcode_handler(
        args.barcode_filename,
        cell_tag=args.cell_tag,
    )

    parse_read = make_parse_read(
        umi_tag=args.umi_tag,
    )

    genotypes = make_genotypes(
        vcf_filename=args.vcf_filename,
        donor_ids=donor_ids,
        prior_strength=args.vcf_prior_strength,
    )

    output_dir = Path(args.output).parent

    if args.genotype_mode == "fixed":
        (
            final_genotypes,
            posterior_probabilities,
            snps,
        ) = run_fixed(
            args,
            genotypes,
            barcode_handler,
            parse_read,
        )

    elif args.genotype_mode == "learn":
        (
            final_genotypes,
            posterior_probabilities,
            snps,
        ) = run_learn(
            args,
            genotypes,
            barcode_handler,
            parse_read,
        )

    elif args.genotype_mode == "refine":
        final_genotypes, posterior_probabilities, snps = run_refine(
            args,
            genotypes,
            barcode_handler,
            parse_read,
        )   

    elif args.genotype_mode == "learn-new-snps":
        (
            final_genotypes,
            posterior_probabilities,
            snps,
        ) = run_learn_new_snps(
            args,
            genotypes,
            barcode_handler,
            parse_read,
            output_dir,
        )

    else:
        raise ValueError(
            f"Unsupported genotype mode: {args.genotype_mode}"
        )

    droplet_type = make_droplet_table(
        posterior_probabilities,
        donor_ids=donor_ids,
        min_prob=args.min_prob,
    )

    snp_summary, fig = summarize_counted_snps(
        snps
    )

    return (
        droplet_type,
        snp_summary,
        fig,
        final_genotypes,
        posterior_probabilities,
    )


def main():
    """Run Demuxalot from the command line."""
    args = parse_args()

    ensure_parent_directory(args.output)
    ensure_bam_index(args.bam_filename)

    (
        droplet_type,
        snp_summary,
        fig,
        final_genotypes,
        posterior_probabilities,
    ) = run_demux(args)

    output_dir = Path(args.output).parent

    droplet_type.reset_index().to_csv(
        args.output,
        sep="\t",
        index=False,
    )

    snp_summary.reset_index().to_csv(
        output_dir / "snp_summary.tsv",
        sep="\t",
        index=False,
    )

    fig.savefig(
        output_dir / "snp_summary.pdf",
        bbox_inches="tight",
    )

    from matplotlib import pyplot as plt

    plt.close(fig)

    if args.posterior_output is not None:
        ensure_parent_directory(
            args.posterior_output
        )

        posterior_probabilities.to_csv(
            args.posterior_output,
            sep="\t",
            index=True,
        )

    if args.learned_genotypes_output is not None:
        ensure_parent_directory(
            args.learned_genotypes_output
        )

        final_genotypes.save_betas(
            args.learned_genotypes_output
        )


if __name__ == "__main__":
    main()
