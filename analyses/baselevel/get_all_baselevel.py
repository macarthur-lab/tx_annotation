from tx_annotation import *

mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_ht_path, "ht")


all_baselevel = get_baselevel_expression_for_genes(mt, gtex, get_proportions=True)

print("writing ht")
all_baselevel = all_baselevel.checkpoint("gs://gnomad-public/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.ht", overwrite=True)
