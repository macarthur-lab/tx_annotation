from tx_annotation import *

manual_curation_file_path = "gs://gnomad-public/papers/2019-tx-annotation/data/gnomad_HIGene_curation.joint.021119.csv"

# Read in and VEP file for tx_annotation
manual_curation = hl.import_table(manual_curation_file_path, delimiter = ",")
manual_curation = manual_curation.annotate(variant_id = manual_curation.variant_id.replace("-",":"))
manual_curation = manual_curation.annotate(** hl.parse_variant(manual_curation.variant_id))
manual_curation = hl.MatrixTable.from_rows_table(manual_curation)
manual_curation = manual_curation.key_rows_by(manual_curation.locus, manual_curation.alleles)
manual_curation = hl.vep(manual_curation, vep_config)

# Tx annotate, reading in the gnomAD mt as a filler so I don't have to write and re-read in manual_curation
filler, gtex = read_tx_annotation_tables(gnomad_release_mt_path,  gtex_v7_tx_summary_mt_path, "ht")
tx_table = tx_annotate_mt(manual_curation, gtex, "proportion", filter_to_csqs=all_coding_csqs)

tx_table.rows().export("gs://gnomad-public/papers/2019-tx-annotation/results/gnomad.HIgene.curation.joint.tx_annotated.021219.tsv.bgz")

