rule sc_match:
    input:

    params:

    output:
        
    container:
        'docker://biagii/scmatch:latest'
    shell:
        'python scMatch.py --refDS refDB/FANTOM5 --dFormat 10x --testDS '


rule to_terms:


rule vis_annos:
