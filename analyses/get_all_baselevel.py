from tx_annotation import *

mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_ht_path, "ht")


all_baselevel = get_baselevel_expression_for_genes(mt, gtex, get_proportions=True)

print("writing ht")
all_baselevel = all_baselevel.checkpoint("gs://gnomad-berylc/tx-annotation/hail2/browser_integration/all.baselevel.021520.ht", overwrite=True)
print("Writing tsv") 
all_baselevel.export("gs://gnomad-berylc/tx-annotation/hail2/browser_integration/all.baselevel.021520.tsv.bgz")
