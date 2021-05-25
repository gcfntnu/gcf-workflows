#!/usr/bin/env python
import sys
import os
import subprocess
import yaml
import glob

def demultiplex_16s_its(args):
    old_dir = os.getcwd()
    os.chdir(args.output_dir)

    primers = get_16s_its_primers(args)
    cutadapt_worker(args, primers)

    os.chdir(old_dir)

def get_16s_its_primers(args):
    with open(args.libprepconf,"r") as libprepconf:
        libconf = yaml.load(libprepconfig, Loader=yaml.Loader)
    l_conf = libconf[args.libkit]
    forward_primers = {}
    reverse_primers = {}
    for region, primers in l_conf['db']['primers'].items():
        forward_primers[region] = primers.split('-')[0]
        reverse_primers[region] = primers.split('-')[1]

    return {'forward': forward_primers, 'reverse': reverse_primers}

def cutadapt_worker(args, primers):
    sample = os.path.basename(args.in1).replace("_R1.fastq.gz","")
    forward = primers.get('forward', {})
    reverse = primers.get('reverse', {})
    os.makedirs("rev_comp_log", exist_ok=True)

    if os.path.exists(os.path.join(args.log_dir, "{}_qiaseq_demultiplex.log".format(sample))):
        cmd = "rm {}".format(os.path.join(args.log_dir, "{}_qiaseq_demultiplex.log".format(sample)))
        subprocess.check_call(cmd, shell=True)

    rev_comp_logs = glob.glob("rev_comp_log/{}_*_rev_comp.log".format(sample))
    if rev_comp_logs:
        cmd = "rm {}".format(" ".join(rev_comp_logs))

    fname = args.in1

    regions = ['unknown']
    for region, primer in forward.items():
        if os.path.exists("{}_unknown_R1.fastq".format(sample)):
            fname = "{}_unknown_R1.fastq".format(sample)
            cp_cmd = "mv -f {} input_{}".format(fname, fname)
            subprocess.check_call(cp_cmd, shell=True)
            cp_cmd = "mv -f {} input_{}".format(fname.replace("R1.fastq","R2.fastq"), fname.replace("R1.fastq","R2.fastq"))
            subprocess.check_call(cp_cmd, shell=True)
            sed_cmd = "sed -i -e s/:region=no_adapter//g input_{}".format(fname)
            subprocess.check_call(sed_cmd, shell=True)
            sed_cmd = "sed -i -e s/:region=no_adapter//g input_{}".format(fname.replace("R1.fastq","R2.fastq"))
            subprocess.check_call(sed_cmd, shell=True)

        overlap_fwd = len(primer.replace("N",""))
        overlap_rev = len(reverse[region].replace("N",""))

        cmd = "cutadapt -g \"{region}={fwd_primer};min_overlap={overlap_fwd};max_error_rate=0.25\" -G \"{region}={rev_primer};min_overlap={overlap_rev};max_error_rate=0.3\" --pair-adapters --untrimmed-output {unknown_r1} --untrimmed-paired-output {unknown_r2} --suffix ':region={{name}}' -o {sample}_{{name}}_R1.fastq -p {sample}_{{name}}_R2.fastq --minimum-length 20 {r1} {r2} >> log/{sample}_qiaseq_demultiplex.log".format(
            sample = sample,
            unknown_r1 = "{}_unknown_R1.fastq".format(sample),
            unknown_r2 = "{}_unknown_R2.fastq".format(sample),
            region = region,
            fwd_primer = primer,
            rev_primer = reverse[region],
            overlap_fwd = overlap_fwd,
            overlap_rev = overlap_rev,
            r1 = ("input_" + fname) if "unknown_R1.fastq" in fname else fname,
            r2 = ("input_" + fname.replace("R1.fastq", "R2.fastq")) if "unknown_R1.fastq" in fname else fname.replace("R1.fastq", "R2.fastq"),
            )
        subprocess.check_call(cmd, shell=True)

        #SECOND PASS OF CUtADAPT TO REMOVE REV COMP PRIMERS IN SHORT AMPLICONS
        if os.path.exists(f"{sample}_{region}_R1.fastq"):
            rev_comp_r = str(Seq(reverse[region]).reverse_complement())
            rev_comp_f = str(Seq(primer).reverse_complement())
            cmd = "cutadapt -a {rev_comp_r} -A {rev_comp_f} --pair-adapters -o tmp_{sample}_{region}_R1.fastq -p tmp_{sample}_{region}_R2.fastq --minimum-length 20 {sample}_{region}_R1.fastq {sample}_{region}_R2.fastq > rev_comp_log/{sample}_{region}_rev_comp.log".format(
                sample = sample,
                region = region,
                rev_comp_r = rev_comp_r,
                rev_comp_f = rev_comp_f,
                )
            subprocess.check_call(cmd, shell=True)

            cmd = f"mv tmp_{sample}_{region}_R1.fastq {sample}_{region}_R1.fastq"
            subprocess.check_call(cmd, shell=True)
            cmd = f"mv tmp_{sample}_{region}_R2.fastq {sample}_{region}_R2.fastq"
            subprocess.check_call(cmd, shell=True)


        regions.append(region)

    if glob.glob("input_{}_*fastq".format(sample)):
        rm_cmd = "rm input_{}_*fastq".format(sample)
        subprocess.check_call(rm_cmd, shell=True)

    R1_comb = ["{}_{}_R1.fastq".format(sample, r) for r in regions]
    R1 = [r1 for r1 in R1_comb if os.path.exists(r1)]
    r1 = " ".join(R1)
    r2 = r1.replace("R1.fastq", "R2.fastq")

    #cat and compress
    cmd = "cat {r1} | pigz -6 -p 2 > {sample}_R1.fastq.gz".format(r1 = r1, sample = sample)
    subprocess.check_call(cmd, shell=True)
    cmd = "cat {r2} | pigz -6 -p 2 > {sample}_R2.fastq.gz".format(r2 = r2, sample = sample)
    subprocess.check_call(cmd, shell=True)

    #unlink r1 and r2 from above
    cmd = "rm -f {}".format(r1)
    subprocess.check_call(cmd, shell=True)
    cmd = "rm -f {}".format(r2)
    subprocess.check_call(cmd, shell=True)
    print("finished sample: {}".format(sample))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--libkit",  help="Library preparation kit name.")
    parser.add_argument("--libprepconfig",  help="Path to libprepconfig.")
    parser.add_argument("--in1",  help="Input R1.")
    parser.add_argument("--in2",  help="Input R2.")
    parser.add_argument("--output-dir",  help="Path to output directory.")
    parser.add_argument("--log-dir",  help="Path to log directory.")

    args = parser.parse_args()

    demultiplex_16s_its(args)
