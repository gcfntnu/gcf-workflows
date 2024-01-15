#!/usr/bin/env python
"""demuxalot pipeline
https://github.com/herophilus/demuxalot/blob/master/examples/2-with-detection-of-new-SNPs.ipynb
"""
import os
import pathlib

import argparse
import pandas as pd
import numpy as np
import pysam

from demuxalot import Demultiplexer, BarcodeHandler, ProbabilisticGenotypes, count_snps

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-b", "--barcode-filename", required=True,
                    help="single cell barcodes file")
parser.add_argument("-i", "--bam-filename", required=True,
                    help="single cell bam file")
parser.add_argument("-v", "--vcf-filename", required=True,
                    help="Reference VCF")
parser.add_argument("-g", "--donor-ids", required=False, default=None,
                    help="comma separated list of donor_ids, or donor_id file. Will fetch donor ids from vcf file as default and will use all samples present in vcf file")
parser.add_argument("-o", "--output", required=True,
                    help="droplet type calls for each barcode")
parser.add_argument("-t", "--threads", default=-1, type=int,
                    help="how many threads to run in parallel (joblib, -1 = all available)")
parser.add_argument("--expected-doublet-rate", default=0.15,
                    help="expected doublet rate.setting to zero skips computation of doublets and helpful when number of clones is too large")


def check_donor_ids_arg(args):
    """validate donor-ids argument
    
    id donor-ids is present -> sort names
    if donor-ids is None > fetch names from vcf. This will use *all* vcf genotype names 
    """
    if args.donor_ids is None:
        with open(args.vcf_filename) as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    header = line.strip().split('\t')
                    genotype_names = list(set(header).difference(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']))
    elif os.path.exists(args.donor_ids):
        # this is a donor_list file (one donor_id per line)
        with open(args.donor_ids) as fh:
            geontype_names = fh.read().splitlines()
    elif "," in args.donor_ids:
        genotype_names = args.donor_ids.split(",")
    else:
        raise ValueError
    
    genotype_names.sort()
    return genotype_names

# copied from demuxalot's utils.py. modified to return fig 
def summarize_counted_SNPs(snp_counts):
    """
    helper function to show number of calls/transcripts available for each barcode
    """
    from collections import Counter
    records = []
    barcode2number_of_calls = Counter()
    barcode2number_of_transcripts = Counter()

    for chromosome, calls in snp_counts.items():
        records.append(dict(
            chromosome=chromosome,
            n_molecules=calls.n_molecules,
            n_snp_calls=calls.n_snp_calls,
        ))

        barcode2number_of_transcripts.update(Counter(calls.molecules["compressed_cb"]))
        barcodes = calls.molecules["compressed_cb"][calls.snp_calls["molecule_index"]]
        barcode2number_of_calls.update(Counter(barcodes))

    from matplotlib import pyplot as plt
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=[12, 5])

    def truncate_at_perc(x, percentile=99.5):
        x = np.asarray(list(x))
        return x.clip(0, np.percentile(x, percentile))

    ax1.hist(truncate_at_perc(barcode2number_of_calls.values()),
             histtype="step", bins=20,
             )
    ax1.set_ylabel("barcodes")
    ax1.set_xlabel("SNP calls per droplet")

    ax2.hist(truncate_at_perc(barcode2number_of_transcripts.values()),
             histtype="step", bins=20,
             )
    ax2.set_ylabel("number of barcodes")
    ax2.set_xlabel("transcripts per droplet")
    #fig.show()

    df = pd.DataFrame(records).sort_values("chromosome").set_index("chromosome")
    
    return df, fig



def run_demux(args):
    handler = BarcodeHandler.from_file(args.barcode_filename)
    genotypes = ProbabilisticGenotypes(genotype_names=args.donor_ids)
    genotypes.add_vcf(args.vcf_filename, prior_strength=100.0)

    snps = count_snps(bamfile_location = args.bam_filename,
                      chromosome2positions=genotypes.get_chromosome2positions(),
                      barcode_handler=handler,
                      joblib_n_jobs=args.threads
                      )

    
    learnt_genotypes, barcode_probs = Demultiplexer.learn_genotypes(snps,
                                                                    genotypes=genotypes,
                                                                    barcode_handler=handler,
                                                                    doublet_prior=args.expected_doublet_rate,
                                                                    )
    
    likelihoods, posterior_probabilities = Demultiplexer.predict_posteriors(snps,
                                                                            genotypes=learnt_genotypes,
                                                                            barcode_handler=handler,
                                                                            doublet_prior=args.expected_doublet_rate,
                                                                            )
    barcode2donor = posterior_probabilities[posterior_probabilities.max(axis=1).gt(0.9)].idxmax(axis=1)
    singlet_pp = posterior_probabilities.loc[:,~posterior_probabilities.columns.str.match("^[A-Za-z0-9_]+\+[A-Za-z0-9_]+$")]
    best_singlet = singlet_pp.idxmax(axis=1)
    droplet_type = pd.DataFrame(index=posterior_probabilities.index)
    droplet_type.index.name = "Barcode"
    droplet_type["doublet_type"] = "unassigned"
    droplet_type.loc[barcode2donor.index] = "singlet"
    doublet_barcodes = barcode2donor.index[barcode2donor.str.match("^[A-Za-z0-9_]+\+[A-Za-z0-9_]+$")]
    droplet_type.loc[doublet_barcodes] = "doublet"
    donor_id = droplet_type.doublet_type.copy()
    singlet_barcodes = droplet_type.index[droplet_type.doublet_type=="singlet"]
    donor_id[singlet_barcodes] = barcode2donor[singlet_barcodes] 
    droplet_type["donor_id"] = donor_id
    droplet_type["best_singlet"] = best_singlet
    snp_summary, fig = summarize_counted_SNPs(snps)
    
    return droplet_type, snp_summary, fig


if __name__ == "__main__":
    args = parser.parse_args()
    args.donor_ids = check_donor_ids_arg(args)
    
    if not pathlib.Path(args.bam_filename + ".bai").exists():
        print("indexing bam file: {} ... \n".format(args.bam_filename))
        pysam.index(bamfile_location)

    droplet_type, snp_summary, fig = run_demux(args)

    droplet_type.reset_index().to_csv(args.output, sep="\t", index=False)
    output_dir = os.path.dirname(args.output)
    snp_summary.reset_index().to_csv(os.path.join(output_dir, "snp_summary.tsv"), sep="\t", index=False)
    fig.savefig(os.path.join(output_dir, "snp_summary.pdf"))
    
