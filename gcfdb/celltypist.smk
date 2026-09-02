rule celltypist_brain_transcriptome_singlecell_atlas_model:
    params:
        proxy = config.get('proxy', {}).get('wget', ''),
        url = 'https://zenodo.org/records/10939707/files/BTS_atlas_CellTypist_model.pkl',
    output:
        model = join(EXT_DIR, 'celltypist', 'data', 'models', 'BTS_atlas_CellTypist_model.pkl')
    threads: 
        24
    log:
        join(EXT_DIR, 'celltypist', 'logs', 'BTS_atlas_CellTypist_model.log')
    shell:
        """
        wget {params.proxy} -O- {params.url} > {output.model}
        echo "An integrative single-cell atlas for exploring the cellular and temporal specificity of genes related to neurological disorders during human brain development. Exp Mol Med 56, 2271–2282 (2024),NA,{params.url},`date -I`" > {log}
        """    


rule celltypist_human_immune_health_atlas_models_L1:
    params:
        proxy = config.get('proxy', {}).get('wget', ''),
        model_url = join(config.get('winecellar', {}).get('url', ''), 'human_immune_health_atlas', 'ref_pbmc_clean_celltypist_model_AIFI_L1_2024-04-18.pkl')
    output:
        model = join(EXT_DIR, 'celltypist', 'data', 'models', 'ref_pbmc_clean_celltypist_model_AIFI_L1_2024-04-18.pkl')
    threads: 
        24
    log:
        join(EXT_DIR, 'celltypist', 'logs', 'ref_pbmc_clean_celltypist_model_AIFI_L1.log')
    shell:
        """
        wget {params.proxy} -O- {params.model_url} > {output.model} 
        echo "Immune Health Atlas 10x 3' CellTypist Model L1,release-2024-04-18,{params.model_url},`date -I`" > {log}
        """

