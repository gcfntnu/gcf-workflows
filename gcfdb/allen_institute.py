
rule mapmycells_prebuild_precomputed_mouse:
    params:
        proxy = config.get('proxy', {}).get('wget', ''),
        url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WMB-10X/20240831/precomputed_stats_ABC_revision_230821.h5',
        release = '20240831'
    output:
        pre_stats = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', 'mus_musculus', 'precomputed_stats.h5')
    threads: 
        24
    log:
        join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells' , 'mus_musculus', 'logs', 'precomputed_stats.log')
    shell:
        """
        wget {params.proxy} -O- {params.url} > {output.pre_stats}
        echo "Allen-Brain-Cell-Atlas,release-{params.release},{params.url},`date -I`" > {log}
        """


rule mapmycells_prebuild_markers_mouse:
    params:
        proxy = config.get('proxy', {}).get('wget', ''),
        url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WMB-10X/20240831/mouse_markers_230821.json',
        release = '2024083'
    output:
        markers = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', 'mus_musculus', 'markers.json')
    threads: 
        24
    log:
        join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells' , 'mus_musculus', 'logs', 'markers.log')
    shell:
        """
        wget {params.proxy} -O- {params.url} > {output.markers}
        echo "Allen-Brain-Cell-Atlas,release-{params.release},{params.url},`date -I`" > {log}
        """


rule mapmycells_prebuild_markers_human:
    params:
        proxy = config.get('proxy', {}).get('wget', ''),
        url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WHB-10Xv3/20240831/query_markers.n10.20240221800.json',
        release = '2024083'
    output:
        markers = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', 'homo_sapiens', 'markers.json')
    threads: 
        24
    log:
        join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells' , 'homo_sapiens', 'logs', 'markers.log')
    shell:
        """
        wget {params.proxy} -O- {params.url} > {output.markers}
        echo "Allen-Brain-Cell-Atlas,release-{params.release},{params.url},`date -I`" > {log}
        """

rule mapmycells_prebuild_precomputed_human:
    params:
        proxy = config.get('proxy', {}).get('wget', ''),
        url = 'https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/mapmycells/WHB-10Xv3/20240831/precomputed_stats.siletti.training.h5',
        release = '20240831'
    output:
        pre_stats = join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells', 'homo_sapiens', 'precomputed_stats.h5')
    threads: 
        24
    log:
        join(EXT_DIR, 'allen-brain-cell-atlas', 'mapmycells' , 'homo_sapiens', 'logs', 'precomputed_stats.log')
    shell:
        """
        wget {params.proxy} -O- {params.url} > {output.pre_stats}
        echo "Allen-Brain-Cell-Atlas,release-{params.release},{params.url},`date -I`" > {log}
        """
