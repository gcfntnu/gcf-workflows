#-*- mode:snakemake -*-
"""Converters between common fileformats
"""

        
rule convert_gtf2txg:
    input:
        join('{ref_dir}', 'anno', 'genes.tsv')
    output:
        join('{ref_dir}', 'anno', 'tx2gene.tsv')
    shell:
        'cut -f1,10 {input} > {output}'
        
rule convert_fasta2sequence_dict:
    input:
        join('{ref_dir}', 'fasta', 'genome.fa')
    output:
        join('{ref_dir}', 'fasta', 'genome.dict')
    singularity:
        'docker://' + config['docker']['picard_gatk']
    shell:
        'picard CreateSequenceDictionary REFERENCE={input} OUTPUT={output}'
        
rule convert_gtf2bed12:
    input:
        join('{ref_dir}', 'anno', 'genes.gtf')
    output:
        join('{ref_dir}', 'anno', 'genes.bed12')
    singularity:
        'docker://' + config['docker']['ucsc-scripts']
    shadow:
        'minimal'
    shell:
        """
        gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons {input} dummy.txt
        genePredToBed dummy.txt {output}
        """
        
rule convert_gtf2intervals:
    input:
        join('{ref_dir}', 'anno', 'genes.gtf')
    output:
        join('{ref_dir}', 'anno', 'genes.intervals')
    shell:
        ''
        
rule convert_gtf2refflat:
    input:
        join('{ref_dir}', 'anno', 'genes.gtf')
    output:
         join('{ref_dir}', 'anno', 'genes.refflat.gz')
    singularity:
        'docker://' + config['docker']['ucsc-scripts']
    shadow:
        'minimal'
    shell:
        """
        gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons {input} dummy.txt
        paste <(cut -f 12 dummy.txt) <(cut -f 1-10 dummy.txt) > refFlat.txt
        gzip refFlat.txt
        mv refFlat.txt.gz {output}
        """
        
rule convert_gtf_transcriptome_gffread:
    input:
        gtf = join('{ref_dir}', 'anno', 'genes.gtf'),
        genome = join('{ref_dir}', 'fasta', 'genome.fa')
    output:
        join('{ref_dir}', 'fasta', 'gtf.gffread.transcripts.fa')
    singularity:
        'docker://' + config['docker']['gffread']
    shell:
        'gffread -w {output} -g {input.genome} {input.gtf}'
        
rule convert_gtf_transcriptome_rsem:
    input:
        gtf = join('{ref_dir}', 'anno', 'genes.gtf'),
        genome = join('{ref_dir}', 'fasta', 'genome.fa')
    params:
        base = join('{ref_dir}', 'fasta', 'gtf.rsem.transcripts.fa')
    singularity:
        'docker://' + config['docker']['rsem']
    output:
        join('{ref_dir}', 'fasta', 'gtf.rsem.transcripts.fa')
    shell:
        'rsem-prepare-reference --gtf {input.gtf} {input.genome} {params.base}'

#fixme
rule convert_transcriptome:
    input:
        join('{ref_dir}', 'fasta', 'gtf.gffread.transcripts.fa')
    output:
        join('{ref_dir}', 'fasta', 'transcriptome.fa')
    shell:
        'ln -sr {input} {output}'

rule convert_gtf2transcript_info:
    input:
        gtf = join('{ref_dir}', 'anno', 'genes.gtf'),
        fasta = join('{ref_dir}', 'fasta', 'transcriptome.fa')
    params:
        script = srcdir('scripts/gtf2tsv.py')
    output:
        tsv = join('{ref_dir}', 'anno', 'transcripts.tsv')
    singularity:
        'docker://' + config['docker']['gcf-bio']
    shell:
        'python {params.script} '
        '--gtf {input.gtf} '
        '--fasta {input.fasta} '
        '--output {output.tsv} '
        '--feature transcript '
        '--add-gc-content '
    
rule exon_gtf:
    input:
       join('{ref_dir}', 'anno', 'genes.gtf')
    output:
        join('{ref_dir}', 'anno', 'exons.gtf')
    singularity:
        'docker://' + config['docker']['default']        
    shell:
        """awk '$3=="exon"' {input} > {output}"""

rule exon_fasta:
    input:
        genome = join('{ref_dir}', 'fasta', 'genome.fa'),
        gtf = join('{ref_dir}', 'anno', 'exons.gtf')
    output:
        join('{ref_dir}', 'fasta', 'exons.fa')
    singularity:
        'docker://' + config['docker']['gcf-bio']        
    shell:
        'bedtools getfasta -fi {input.genome} -fo {output} -bed {input.gtf} '    

        
rule convert_gtf2gene_info:
    input:
        gtf = join('{ref_dir}', 'anno', 'genes.gtf'),
        fasta = join('{ref_dir}', 'fasta', 'exons.fa')
    params:
        script = srcdir('scripts/gtf2tsv.py'),
    output:
        tsv = join('{ref_dir}', 'anno', 'genes.tsv')
    singularity:
        'docker://' + config['docker']['gcf-bio']
    shell:
        'python {params.script} '
        '--gtf {input.gtf} '
        '--fasta {input.fasta} '
        '--output {output.tsv} '
        '--feature gene '
        '--add-gc-content '
        
