from tx_annotation import *

context_ht_path = 'gs://gnomad-resources/context/hail-0.2/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht/'
mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_mt_path, "ht")


all_baselevel = get_baselevel_expression_for_genes(mt, gtex, get_proportions=True)

print("Writing tsv") 
all_baselevel.export("gs://gnomad-berylc/tx-annotation/hail2/browser_integration/all.baselevel.021519.tsv.bgz")
print("writing ht") 
all_baselevel.write("gs://gnomad-berylc/tx-annotation/hail2/browser_integration/all.baselevel.021519.ht", overwrite=True)
