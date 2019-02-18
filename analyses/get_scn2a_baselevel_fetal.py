from tx_annotation import *
mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_mt_path, "ht")

gene_baselevel= get_baselevel_expression_for_genes(mt, gtex, gene_list = {'SCN2A'})
gene_baselevel.export("gs://gnomad-public/papers/2019-tx-annotation/results/SCN2A.baselevel.ext.021719.tsv.bgz")

gene_baselevel_prop= get_baselevel_expression_for_genes(mt, gtex, gene_list = {'SCN2A'}, get_proportions=True)
gene_baselevel_prop.export("gs://gnomad-public/papers/2019-tx-annotation/results/SCN2A.baselevel.pext.021719.tsv.bgz")


hbdr_fetal_path = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR.RSEM.sample_specific.tx_medians.021719.mt"
mt, hbdr_fetal = read_tx_annotation_tables(context_ht_path, hbdr_fetal_path, "ht")
gene_baselevel= get_baselevel_expression_for_genes(mt, hbdr_fetal, gene_list = {'SCN2A'})
gene_baselevel.export("gs://gnomad-public/papers/2019-tx-annotation/results/SCN2A.HBDR.fetal.baselevel.ext.02171