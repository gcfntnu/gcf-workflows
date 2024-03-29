#-*- mode:snakemake -*-
from os.path import dirname, abspath

include:
    join(GCFDB_DIR, 'quickmirseq.db')
    
rule quickmirseq_symlink:
    input:
        unpack(get_filtered_fastq)
    params:
        fastq_dir = join(QUANT_INTERIM, 'quickmirseq', 'fastq'),
        samples = join(QUANT_INTERIM, 'quickmirseq', 'allids.txt')
    output:
        fastq = join(QUANT_INTERIM, 'quickmirseq', 'fastq', '{sample}.fastq')
    run:
        import os
        src = abspath(input.R1)
        dst = abspath(output.fastq)
        os.symlink(src, dst)
        #print('symlinking: {} > {}'.format(src, dst))
        with open(params.samples, 'a') as fh:
            name, ext = os.path.splitext(os.path.basename(dst))
            fh.write(name + '\n')
        
            
rule quickmirseq_config:
    input:
        db = join(EXT_DIR, 'quickmirseq', config['organism'], '.done')
    params:
        fastq_dir = rules.quickmirseq_symlink.params.fastq_dir,
        species = config['organism'],
        rna_db = abspath(join(EXT_DIR, 'quickmirseq', ORG)),
        dna_bowtie_index = abspath(join(EXT_DIR, 'quickmirseq', ORG, 'genome')),
        output = join(QUANT_INTERIM, 'quickmirseq')
    threads:
        24
    output:
        join(QUANT_INTERIM, 'quickmirseq', 'quickmirseq.conf')
    shell:
        """
        echo SPECIES={params.species} >> {output}
        echo RNA_BOWTIE_INDEX={params.rna_db} >> {output}
        echo GENOME_BOWTIE_INDEX={params.dna_bowtie_index} >> {output}
        echo FASTQ_DIR={params.fastq_dir} >> {output}
        echo CPU={threads} >> {output}
        echo STRAND=1 >> {output}
        echo FASTQ_SUFFIX=fastq >> {output}
        echo OUTPUT_FOLDER={params.output} >> {output}
        echo RUN_EXTENSION_EVALUATION=yes >> {output}
        echo CUTADAPT_REQUIRED=no >> {output}
        echo REFINE_MISMATACH_READS=yes >> {output}
        echo ZEROCOUNT_THRESHOLD=1 >> {output}
        echo ZEROCOUNT_SAMPLE_THRESHOLD=0.6 >> {output}
        echo AVG_READ_THRESHOLD=2 >> {output}
        echo MIN_MIRNA_LENGTH=16 >> {output}
        echo MAX_MIRNA_LENGTH=42 >> {output}
        echo EXTENSION3=5 >> {output}
        echo EXTENSION5=4 >> {output}
        echo keepTemp=yes >> {output}
        echo UNIQUE_LIBRARY_ANALYSIS=yes >> {output}
        echo OFFSET_ANALYSIS=yes >> {output}
        """    

rule quickmirseq_quant:
    input:
        samples = expand(join(QUANT_INTERIM, 'quickmirseq', 'fastq', '{sample}.fastq'), sample=SAMPLES),
        config = join(QUANT_INTERIM, 'quickmirseq', 'quickmirseq.conf')
    params:
        fastq_dir = rules.quickmirseq_symlink.params.fastq_dir,
        sample_file = join(QUANT_INTERIM, 'quickmirseq', 'allids.txt')
    output:
        mir_table = join(QUANT_INTERIM, 'quickmirseq', 'miR.Counts.csv'),
        isomir_table = join(QUANT_INTERIM, 'quickmirseq', 'isoform.Counts.csv'),
        annotation = join(QUANT_INTERIM, 'quickmirseq', 'readAnnotDistribution.csv'),
        unmapped = join(QUANT_INTERIM, 'quickmirseq', 'unmapped.csv')
    threads:
        8
    container:
        'docker://' + config['docker']['quickmirseq']
    shell:
        'perl  $QuickMIRSeq/QuickMIRSeq.pl  '
        '{params.sample_file} '
        '{input.config} '

rule quickmirseq_report:
    input:
        rules.quickmirseq_quant.output
    params:
        outdir = join(QUANT_INTERIM, 'quickmirseq')
    output:
        join(QUANT_INTERIM, 'quickmirseq', 'miRNA-reads.json')
    container:
        'docker://' + config['docker']['quickmirseq']
    shell:
        """
        cd {params.outdir}
        QuickMIRSeq-report.sh
        """

rule quickmirseq_all:
    input:
        rules.quickmirseq_report.output
