#!/usr/bin/env python
"""QIIME2 q2 classifier
"""

import sys
import argparse
from os.path import join, dirname, basename, splitext

from qiime2 import Artifact
from qiime2.plugins.feature_classifier import methods


def load_data(rep_set, taxa):
    q2_repset = Artifact.import_data('FeatureData[Sequence]', rep_set)
    q2_taxa = Artifact.import_data('FeatureData[Taxonomy]', taxa, view_type='HeaderlessTSVTaxonomyFormat')
    return q2_repset, q2_taxa

def extract_region_fasta(q2_repset, fwd='CCTACGGGNGGCWGCAG', rev='GACTACHVGGGTATCTAATCC', min_len=None, n_jobs=1):
    q2_refseq_filtered = methods.extract_reads(sequences=q2_repset, f_primer=fwd, r_primer=rev, min_length=min_len, n_jobs=n_jobs)
    return q2_refseq_filtered

def train_classifier(q2_refseq, q2_taxa):
    # train the classifier
    clf = methods.fit_classifier_naive_bayes(reference_reads=q2_refseq, reference_taxonomy=q2_taxa)
    return clf



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', help='rep set sequences', type=argparse.FileType('r'))
    parser.add_argument('--taxa', help='taxonomy files', type=argparse.FileType('r'))
    parser.add_argument('--output', help="classifier output file (.qza)", default='classifier.qza')
    parser.add_argument('--f-primer', help="filer reference fasta wrt forward primer compatibility")
    parser.add_argument('--r-primer', help="filer reference fasta wrt forward reverse compatibility")
    parser.add_argument('--min-len', help="Minimum amplicon length", type=int, default=30)
    parser.add_argument('--threads', help="Number of threads to use", type=int, default=1)
    args = parser.parse_args(sys.argv[1:])

    q2_repset, q2_taxa = load_data(args.fasta.name, args.taxa.name)
    q2_refseq_filtered = extract_region_fasta(q2_repset, fwd=args.f_primer, rev=args.r_primer, min_len=args.min_len, n_jobs=args.threads)
    classifier = train_classifier(q2_refseq_filtered.reads, q2_taxa)
    classifier.classifier.save(args.output)
    
    classifier.classifier.export_data(join(dirname(args.output), 'export', splitext(basename(args.output))[0] ))
