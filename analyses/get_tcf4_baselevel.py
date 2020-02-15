'''Code to create TCF4 plot in Fig2B'''

from tx_annotation import *

context_ht_path = 'gs://gnomad-resources/context/hail-0.2/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht'

'''
Fig2B, SuppFig3
mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_mt_path, "ht")

gene_baselevel= get_baselevel_expression_for_genes(mt, gtex, gene_list = {'TCF4'})
gene_baselevel.export("gs://gnomad-public/papers/2019-tx-annotation/results/TCF4.baselevel.ext.021319.tsv.bgz")

gene_baselevel_prop= get_baselevel_expression_for_genes(mt, gtex, gene_list = {'TCF4'}, get_proportions=True)
gene_baselevel_prop.export("gs://gnomad-public/papers/2019-tx-annotation/results/TCF4.baselevel.pext.021319.tsv.bgz")

'''
'''
salmon_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/salmon.GTEx.v7.brain.cortex.tx_medians.021519.mt"
mt, salmon = read_tx_annotation_tables(context_ht_path, salmon_brain_cortex_out_path, "ht")
salmon_baselevel_prop= get_baselevel_expression_for_genes(mt, salmon, gene_list = {'TCF4'}, get_proportions=True)
salmon_baselevel_prop.export("gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/TCF4/salmon.TCF4.baselevel.pext.021519.tsv.bgz")

rsem_brain_cortex_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.GTEx.v7.brain.cortex.tx_medians.021519.mt"
mt, rsem = read_tx_annotation_tables(context_ht_path, rsem_brain_cortex_path, "ht")
rsem_baselevel_prop= get_baselevel_expression_for_genes(mt, rsem, gene_list = {'TCF4'}, get_proportions=True)
rsem_baselevel_prop.export("gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/TCF4/rsem.TCF4.baselevel.pext.021519.tsv.bgz")
'''


mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_mt_path, "ht")

gene_baselevel= get_baselevel_expression_for_genes(mt, gtex, gene_list = {'SCN2A'})
gene_baselevel.export("gs://gnomad-public/papers/2019-tx-annotation/results/SCN2A.baselevel.ext.021719.tsv.bgz")

gene_baselevel_prop= get_baselevel_expression_for_genes(mt, gtex, gene_list = {'SCN2A'}, get_proportions=True)
gene_baselevel_prop.export("gs://gnomad-public/papers/2019-tx-annotation/results/SCN2A.baselevel.pext.021719.tsv.bgz")


hbdr_fetal_path = "gs://gnomad-berylc/tx-annotation/hail2/fetal_rnaseq/manuscript/HBDR.RSEM.sample_specific.tx_medians.021719.mt"
mt, hbdr_fetal = read_tx_annotation_tables(context_ht_path, hbdr_fetal_path, "ht")
gene_baselevel= get_baselevel_expression_for_genes(mt, hbdr_fetal, gene_list = {'SCN2A'})
gene_baselevet("gs://gnomad-public/papers/2019-tx-annotation/results/SCN2A.HBDR.fetal.baselevel.ext.021719.tsv.bgz")l.expor