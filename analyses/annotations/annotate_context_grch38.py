from tx_annotation import *

v8_isoform_expression = "gs://gnomad-public/papers/2019-tx-annotation/data/GRCH38/" \
                        "reheadered.GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"

gtex_v8_tx_summary_ht_path = "gs://gnomad-public/papers/2019-tx-annotation/data/GRCH38" \
                             "/GTEx.V8.tx_medians.021620.ht"
gtex_v8_gene_maximums_ht_path = "gs://gnomad-public/papers/2019-tx-annotation/data/GRCH38/" \
                                "GTEx.v8.gene_expression_per_gene_per_tissue.021620.ht"

context_hg38_annotated = "gs://gnomad-public/papers/2019-tx-annotation/data/GRCH38/all.possible.snvs.tx_annotated.hg38.021620.ht"
context_hg38_max_per_gene= "gs://gnomad-public/papers/2019-tx-annotation/data/GRCH38/all.genes.max.pext.hg38.021620.tsv.bgz"

# make isoform summary file
get_gtex_summary(v8_isoform_expression, gtex_v8_tx_summary_ht_path)

# make maximum expression file summary file
get_gene_expression(gtex_v8_tx_summary_ht_path, gtex_v8_gene_maximums_ht_path)

# Annotate context (ie. all possible bases)
mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v8_tx_summary_ht_path, 'ht')

mt_annotated = tx_annotate_mt(mt, gtex, "proportion",
                              tissues_to_filter=None,
                              gene_maximums_ht_path=gtex_v8_gene_maximums_ht_path,
                              filter_to_csqs=all_coding_csqs)

mt_annotated.write(context_hg38_annotated, overwrite=True)

identify_maximum_pext_per_gene(context_hg38_annotated, context_hg38_max_per_gene)