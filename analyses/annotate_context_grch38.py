from tx_annotation import *


hg38_tx_summary_mt_path = "gs://gnomad-berylc/tx-annotation/hail2/GTEx.V8.tx_medians.020920.mt"
hg38_gene_maximums_kt_path = "gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.v8.gene_expression_per_gene_per_tissue.020920.kt"

context_hg38_annotated = "gs://gnomad-public/papers/2019-tx-annotation/data/all.possible.snvs.tx_annotated.HBDR.020920.ht"
context_hg38_max_per_gene = "gs://gnomad-public/papers/2019-tx-annotation/data/all.genes.max.pext.HBDR.020920.tsv.bgz"


mt, hg38 = read_tx_annotation_tables(context_ht_path,hg38_tx_summary_mt_path, 'ht')

print("Starting tx annotation on")

mt_annotated = tx_annotate_mt(mt, hg38, "proportion",
                              tissues_to_filter=None,
                              gene_maximums_kt_path = hg38_gene_maximums_kt_path,
                              filter_to_csqs=all_coding_csqs)

print("Finished annotation ")
#
mt_annotated.write(context_hg38_annotated, overwrite=True)
#
print("Wrote RSEM Brain Cortex HT, now going to identify maximums")
#
identify_maximum_pext_per_gene(context_hg38_annotated, context_hg38_max_per_gene)