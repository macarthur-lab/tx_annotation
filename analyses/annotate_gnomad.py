from tx_annotation import *
mt, gtex = read_tx_annotation_tables(gnomad_release_mt_path, gtex_v7_tx_summary_ht_path, 'ht')

mt_annotated = tx_annotate_mt(mt, gtex,"proportion",
                              filter_to_csqs=all_coding_csqs)
mt_annotated.write("gs://gnomad-public/papers/2019-tx-annotation/data/gnomad_release_annotated/gnomad.exomes.r2.1.1.sites.tx_annotated.021520.ht")
