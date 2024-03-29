#-*- mode:snakemake -*-

ASSEMBLY_INTERIM = join(INTERIM_DIR, 'rnaseq', 'assembly')

include:
    join(GCFDB_DIR, 'reference_db.db')
include:
   join(GCFDB_DIR,  'spikein.db')
include:
    join(GCFDB_DIR, 'contaminants.db')
include:
    join(GCFDB_DIR, 'rrna.db')

BASE_REF_DIR = config['base_ref_dir']

if config['assembly']['transcriptome']['assembler'] != 'skip' or config['assembly']['genome']['assembler'] != 'skip':
    REF_DIR = ASSEMBLY_INTERIM
    include:
        'assembly/align.smk'
    include:
        'assembly/transcriptome.smk'
    include:
        'assembly/genome.smk'


    def get_tx_assembly(*args):
        merger = config['assembly']['transcriptome'].get('merger', 'skip')
        assembler = config['assembly']['transcriptome'].get('assembler', 'skip')
        if merger == 'skip' or merger == 'skip':
            gtf = join(BASE_REF_DIR, 'anno', 'genes.gtf')
            fasta = join(BASE_REF_DIR, 'fasta', 'gtf.rsem.transcripts.fa')
            return {'gtf': gtf, 'fasta': fasta}

        if assembler == 'scallop' and merger == 'taco':
            gtf = join(ASSEMBLY_INTERIM, 'taco', 'scallop', 'assembly.refcomp.gtf')
            fasta = join(ASSEMBLY_INTERIM, 'taco', 'scallop', 'transcriptome.fa'),
        elif assembler == 'stringtie' and merger == 'taco':
            gtf = join(ASSEMBLY_INTERIM, 'taco', 'stringtie', 'assembly.refcomp.gtf'),
            fasta = join(ASSEMBLY_INTERIM, 'taco', 'stringtie', 'transcriptome.fa')
        elif assembler == 'stringtie' and merger == 'stringtie_merge':
            gtf = join(ASSEMBLY_INTERIM, 'stringtie_merge', 'stringtie', 'assembly.refcomp.gtf'),
            fasta = join(ASSEMBLY_INTERIM, 'stringtie_merge', 'stringtie', 'transcriptome.fa')
        elif assembler == 'scallop' and merger == 'stringtie_merge':
            gtf = join(ASSEMBLY_INTERIM, 'stringtie_merge', 'scallop', 'assembly.refcomp.gtf'),
            fasta = join(ASSEMBLY_INTERIM, 'stringtie_merge', 'scallop', 'transcriptome.fa')
        elif assembler == 'psiclass' and merger == 'psiclass':
            gtf = join(ASSEMBLY_INTERIM, 'psiclass', 'anno', 'assembly.refcomp.gtf'),
            fasta = join(ASSEMBLY_INTERIM, 'psiclass', 'fasta', 'transcriptome.fa')
        else:
            print(config['assembly'])
            raise ValueError("Wrong assembler config")
        return {'gtf': gtf, 'fasta': fasta}

    rule assembled_transcriptome:
        input:
            unpack(get_tx_assembly)
        output:
            gtf = join(REF_DIR, 'anno', 'genes.gtf'),
            fasta = join(REF_DIR, 'fasta', 'transcriptome.fa')
        shell:
            """
            ln -sr {input.gtf} {output.gtf} 
            ln -sr {input.fasta} {output.fasta} 
            """

    def get_genome_assembly(*args):
        assembly_method = config['assembly']['genome'].get('method', 'skip')
        if not assembly_method == 'skip':
            raise NotImplementedError

        return {'genome': join(BASE_REF_DIR, 'fasta', 'genome.fa'),
                'genome_index': join(BASE_REF_DIR, 'fasta', 'genome.fa.fai'),
                'index': directory(join(BASE_REF_DIR, 'index'))}

    rule assembled_genome:
        input:
            unpack(get_genome_assembly)
        output:
            genome = join(REF_DIR, 'fasta', 'genome.fa'),
            genome_index = join(REF_DIR, 'fasta', 'genome.fa.fai')
        params:
            fasta_dir = join(REF_DIR, 'fasta'),
            index_dir = join(REF_DIR, 'index')
        shell:
            """
            mkdir -p {params.fasta_dir}
            mkdir -p {params.index_dir}
            ln -srf {input.genome} {output.genome}
            ln -srf {input.genome_index} {output.genome_index}
            ln -srf {input.index}/genome {params.index_dir}/genome
            """
