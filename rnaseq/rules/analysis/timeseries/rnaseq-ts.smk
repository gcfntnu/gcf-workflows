TS_INTERIM = join(ANALYSIS_INTERIM, 'timeseries')


rule ts_inspect_quant:
    input:
        bam = expand(rules.star_align.output.bam, sample=SAMPLES),
        txdb = 
    script:
        srcdir('scripts/inspect_minus_ts.R')
    params:
        orientation = config['read_geometry'],
        
    
