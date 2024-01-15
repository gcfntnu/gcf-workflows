#! /usr/bin/env python
"""sort vcf header lines to match a header lines of a bam file

"""
import sys
import os
import re
import subprocess
from pathlib import Path, PurePath

import argparse
import pandas as pd


def get_bam_chromosomes(bam_fn):
    result = subprocess.run(['samtools', 'view', '-H', bam_fn], stdout=subprocess.PIPE) 
    if result.returncode:
        print(f'Command {result.cmd} failed with error {result.returncode}')
    else:
        out = []
        for line in result.stdout.decode().splitlines():
            if line.startswith('@SQ'):
                ind, sn, ln = line.split('\t')
                out.append([sn.split(':')[-1], ln.split(':')[-1]])
        return pd.DataFrame(out, columns=['SN', 'LN'])
    
def extract_vcf_header(fh):
    contigs_found = False
    pre_lines, contig_lines, post_lines = [],[],[]
    
    for line in fh:
        if not contigs_found and not line.startswith('##contig'):
            pre_lines.append(line)
        elif line.startswith('##contig'):
            contigs_found = True
            contig_lines.append(line)
        else:
            if not line.startswith('#CHROM') and contigs_found:
                post_lines.append(line)
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            break
    return pre_lines, contig_lines, post_lines, fh, header

def get_vcf_chromsomes(contig_lines):
    ll = []
    for line in contig_lines:
        els = re.findall('##contig=<(.*)>$', line.strip())[0].split(',')
        series = pd.DataFrame.from_dict(dict(x.split("=") for x in els), orient='index').T
        ll.append(series)
    df = pd.concat(ll, axis=0)
    if 'ID' in df.columns:
        if 'length' in df.columns:
            return df[['ID', 'length']]
        return df[['ID']]
    else:
        raise ValueError('failed to find ID and length keywords in vcf file')
    
def to_rows(anno):
    rowdicts = []
    failed = []
    for i, line in enumerate(anno.values):
        lx = (it.split("=") if "=" in it  else [it, None] for it in line.split(";"))
        rowdicts.append({k: v for k, v in lx})
    return pd.DataFrame.from_records(rowdicts).set_index(anno.index)

def read_vcf(f, names=None, chunksize=int(1e4), extended=False, comment=None, header_chr_encoded=False):
    dtypes = {"#CHROM": "category", "FILTER": "category"}
    df_iter = pd.read_csv(f,
                          sep="\t",
                          header=None,
                          names=names,
                          dtype=dtypes,
                          chunksize=chunksize,
                          comment=comment
                          )
    dfs = []
    for df in df_iter:
        if extended:
            cols = df.columns
            extra = to_rows(df.INFO.astype(str))
            df = df.drop("INFO", axis=1)
            extra.set_index(df.index, inplace=True)
            try:
                extra = extra.astype('f')
            except:
                pass
            ndf = pd.concat([df, extra], axis=1, sort=False)
            chrom_encoded = ndf['#CHROM'].iloc[:10].str.startswith('chr').sum() == 10
            if header_chr_encoded and not chrom_encoded:
                ndf['#CHROM'] = 'chr' + ndf['#CHROM'].astype('string')
            elif chrom_encoded and not header_chr_encoded:
                ndf['#CHROM'] =  [i.split('chr')[-1] for i in  ndf['#CHROM']]
            dfs.append(ndf)
        else:
            dfs.append(df)

    df = pd.concat(dfs, sort=False)
    df.index = list(df['#CHROM'] + ':' + df['POS'].astype('string'))

    if extended:
        return df, cols
    return df


def new_contig_lines(chroms_ordered, vcf_uniq_chroms):
    sorted_lines = []
    for name in chroms_ordered:
        line = '##contig=<ID={}>\n'.format(name)
        sorted_lines.append(line)
    for name in vcf_uniq_chroms:
        sorted_lines.append(line)
    return sorted_lines

def write_vcf(output_fn, header_lines, vcf):
    with open(output_fn, 'w') as fh:
        for line in header_lines:
            fh.write(line)
        vcf.to_csv(fh, index=False, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf-file', help='input vcf file')
    parser.add_argument('--bam-file', help='input bam file')
    parser.add_argument('--output-file', default="common_variants_grch38.sorted.vcf.gz", help='output vcf file. if .bgz/.gz ending it will be compressed with bgzip and indexed by tabix')
    parser.add_argument('--write-ids-file', help='write chrom-start-end formatted file (scsplit)', action='store_true')
    parser.add_argument('--filter-eur-af', default=None, help="filter vcf field EUR_AF > filter_eur_af")
    parser.add_argument('--rename-chrom-file', default=None, help="rename chromsome names in vcf file with this mapping file (tab separated, two columns)")
    
    
    args = parser.parse_args()
    args.vcf_file = Path(args.vcf_file)
    args.bam_file = Path(args.bam_file)
    args.output_file = PurePath(args.output_file)

    assert args.vcf_file.suffix == '.vcf'
    assert args.bam_file.suffix == '.bam'
    if args.output_file.name.endswith('.gz') or args.output_file.name.endswith('.bgz'):
        suffix = ''.join(args.output_file.suffixes)
        args.output_file = args.output_file.with_suffix('')
        compress = True
    
    bam_chroms = get_bam_chromosomes(args.bam_file)
    bam_chr_encoded = bam_chroms.SN.str.contains('chr').sum() > 20
    with open(args.vcf_file) as fh:
        vcf_header_pre_lines, vcf_contig_lines, vcf_header_post_lines, fh, cols = extract_vcf_header(fh)
        vcf_chroms = get_vcf_chromsomes(vcf_contig_lines)
        vcf_chr_encoded = vcf_chroms.ID.str.contains('chr').sum() > 20
        if args.filter_eur_af is not None:
            vcf, cols = read_vcf(fh, cols, header_chr_encoded=vcf_chr_encoded, extended=True)
        else:
            vcf = read_vcf(fh, cols, header_chr_encoded=vcf_chr_encoded)

    if not bam_chr_encoded and vcf_chr_encoded:
        new_ids = vcf_chroms.ID.str.extract('chr(.*)')
        vcf_chroms['ID'] = new_ids

    if bam_chr_encoded and not vcf_chr_encoded:
        new_ids = 'chr' + vcf_chroms.ID
        vcf_chroms['ID'] = new_ids

    vcf_ID_in_bam = vcf_chroms['ID'].isin(bam_chroms['SN'])
    keep = vcf_chroms.loc[vcf_ID_in_bam,]

    # list of chromosomes ordered by bam file header
    chroms_ordered = [i for i in bam_chroms['SN'] if i in keep['ID'].values]
    vcf_uniq_chroms = [] 
    if not all(vcf_ID_in_bam): #if we have chromosome names in vcf not present in bam
        vcf_uniq_chroms = vcf_chroms.loc[~vcf_ID_in_bam,].to_list()
        keep_chroms = chroms_ordered
        vcf_chrom_chrcoded = sum(["chr" in i for i  in vcf['#CHROM'].unique()]) > 20
        if bam_chr_encoded and not vcf_chrom_chrcoded:
            keep_chroms = [i.split('chr')[-1] for i in chroms_ordered]
        vcf = vcf[vcf['#CHROM'].isin(keep_chroms)]
        
    if args.filter_eur_af:
        vcf_filtered = vcf[vcf.EUR_AF.astype('d') > args.filter_eur_af]
        with open(args.vcf_file) as fh:
            vcf = read_vcf(fh, cols, header_chr_encoded=vcf_chr_encoded)
        vcf = vcf.loc[vcf_filtered.index,:]

    if bam_chr_encoded and (vcf['#CHROM'].iloc[:10].str.startswith('chr').sum() == 0):
        vcf['#CHROM'] = 'chr' + vcf['#CHROM'].astype('string')
        vcf.index = list(vcf['#CHROM'] + ':' + vcf['POS'].astype('string'))      
    assert all(vcf['#CHROM'].isin(chroms_ordered))
    
    sorted_contig_lines = new_contig_lines(chroms_ordered, vcf_uniq_chroms)
    header_lines = vcf_header_pre_lines
    header_lines.extend(sorted_contig_lines)
    header_lines.extend(vcf_header_post_lines)
    write_vcf(args.output_file, header_lines, vcf)
    
    if args.write_ids_file: #scsplit format
        snp_ids = pd.DataFrame(vcf.index)
        p = Path(args.output_file)
        ids_output = p.with_suffix('.ids')
        snp_ids.to_csv(ids_output, sep="\t", index=False, header=False)

    if compress:
        # write index together with compressed file
        index_suffix = suffix + '.tbi'
        cmd = 'bgzip --threads 4 --stdout --index --index-name {} {} > {}'.format(str(args.output_file.with_suffix(index_suffix)), str(args.output_file), str(args.output_file.with_suffix(suffix)))
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)
        
    if args.rename_chrom_file is not None:
        # bcftools annotate --rename-chrs chr_name_conv.txt in.vcf -Oz -o out.gcf
        raise NotImplementedError
