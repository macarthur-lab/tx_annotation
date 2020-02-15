from tx_annotation import *

rsem_tpms = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.RSEM.tissue_names.tpm.tsv.gz"

get_gtex_summary(rsem_tpms, "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.tissue_names.medians.021520.mt")

get_gene_expression("gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.tissue_names.medians.021520.mt",
                    "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR.gene_expression_per_gene_per_tissue.021520.ht")


# At this point updated hbdr_fetal_tissue_summary_ht_path, hbdr_fetal_tissue_gene_maximums_ht_path in tx_annotation_resources.py
print(hbdr_tx_summary_mt_path)
print(hbdr_gene_maximums_ht_path)

context_hbdr_annotated = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/all.possible.snvs.tx_annotated.HBDR.021520.ht"
context_hbdr_max_per_gene = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/all.genes.max.pext.HBDR.021520.tsv.bgz"

mt, hbdr = read_tx_annotation_tables(context_ht_path, hbdr_tx_summary_mt_path, 'ht')

print("Starting tx annotation on")

mt_annotated = tx_annotate_mt(mt, hbdr, "proportion",
                              tissues_to_filter=None,
                              gene_maximums_ht_path=hbdr_gene_maximums_ht_path,
                              filter_to_csqs=all_coding_csqs)

print("Finished annotation ")

mt_annotated.write(context_hbdr_annotated, overwrite=True)

print("Wrote RSEM Brain Cortex HT, now going to identify maximums")

identify_maximum_pext_per_gene(context_hbdr_annotated, context_hbdr_max_per_gene)