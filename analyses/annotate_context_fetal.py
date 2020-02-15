from tx_annotation import *


hbdr_tx_summary_mt_path = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.tissue_names.medians.mt"
hbdr_gene_maximums_kt_path = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR.gene_expression_per_gene_per_tissue.020920.kt"

context_hbdr_annotated = "gs://gnomad-public/papers/2019-tx-annotation/data/all.possible.snvs.tx_annotated.HBDR.020920.ht"
context_hbdr_max_per_gene = "gs://gnomad-public/papers/2019-tx-annotation/data/all.genes.max.pext.HBDR.020920.tsv.bgz"


mt, hbdr = read_tx_annotation_tables(context_ht_path,hbdr_tx_summary_mt_path, 'ht')

print("Starting tx annotation on")

mt_annotated = tx_annotate_mt(mt, hbdr, "proportion",
                              tissues_to_filter=None,
                              gene_maximums_kt_path = hbdr_gene_maximums_kt_path,
                              filter_to_csqs=all_coding_csqs)

print("Finished annotation ")
#
mt_annotated.write(context_hbdr_annotated, overwrite=True)
#
print("Wrote RSEM Brain Cortex HT, now going to identify maximums")
#
identify_maximum_pext_per_gene(context_hbdr_annotated, context_hbdr_max_per_gene)