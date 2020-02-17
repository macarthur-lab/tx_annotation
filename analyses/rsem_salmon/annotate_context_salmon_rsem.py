'''
In order to stay as consistent as possible in the figure regeneration with removing low best genes. Need to
1 - Annotate context ht path with 151 Brain Cortex samples with RSEM
2 - Annotate context ht path with 151 Brain Cortex samples with Salmon
3 - Get maximums for both
4 - Filter maximums in ClinVar analysis
5 - Filter maximums in de novo analysis
'''

from tx_annotation import *


rsem_brain_cortex_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.GTEx.v7.brain.cortex.tx_medians.021519.mt"
#context_rsem_brain_cortex_annotated = "gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/gene_list_comparisons_salmon_rsem/all.possible.snvs.tx_annotated.rsem.brain.cortex.020420.ht"
#context_rsem_brain_cortex_max_per_gene= "gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/gene_list_comparisons_salmon_rsem/all.genes.max.pext.rsem.brain.cortex.020420.tsv.bgz"


salmon_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/salmon.GTEx.v7.brain.cortex.tx_medians.021519.mt"
#context_salmon_brain_cortex_annotated = "gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/gene_list_comparisons_salmon_rsem/all.possible.snvs.tx_annotated.salmon.brain.cortex.020420.ht"
#context_salmon_brain_cortex_max_per_gene= "gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/gene_list_comparisons_salmon_rsem/all.genes.max.pext.salmon.brain.cortex.020420.tsv.bgz"


# 1 - Annotate context ht path with RSEM Brain Cortex
# mt, rsem = read_tx_annotation_tables(context_ht_path,rsem_brain_cortex_path, 'ht')
#
# print("Starting tx annotation on")
#
# mt_annotated = tx_annotate_mt(mt, rsem, "proportion",
#                               tissues_to_filter=None,
#                               gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path,
#                               filter_to_csqs=all_coding_csqs)
#
# print("Finished annotation ")
#
# mt_annotated.write(context_rsem_brain_cortex_annotated, overwrite=True)
#
# print("Wrote RSEM Brain Cortex HT, now going to identify maximums")
#
# identify_maximum_pext_per_gene(context_rsem_brain_cortex_annotated, context_rsem_brain_cortex_max_per_gene)
#
# print("Done with RSEM Brain - Cortex")


mt, rsem = read_tx_annotation_tables(gnomad_release_mt_path, rsem_brain_cortex_path, 'ht')

mt_annotated = tx_annotate_mt(mt, rsem,"proportion", tissues_to_filter=None, gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path,filter_to_csqs=all_coding_csqs)
mt_annotated.write("gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/data/gnomad.exomes.r2.1.1.sites.tx_annotated.brain.cortex.RSEM.020520.mt",overwrite = True)


# Annotate gnomAD sites MT (I had already done this a year ago but want to make sure it's correct so doing it again)


# 2 - Annotate context ht path with salmon Brain Cortex
# mt, salmon = read_tx_annotation_tables(context_ht_path, salmon_brain_cortex_out_path, 'ht')
#
# print("Starting tx annotation on")
# mt_annotated = tx_annotate_mt(mt, salmon, "proportion",
#                               tissues_to_filter=None,
#                               gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path,
#                               filter_to_csqs=all_coding_csqs)
# print("Finished annotation ")
#
# mt_annotated.write(context_salmon_brain_cortex_annotated, overwrite=True)
#
# print("Wrote Salmon Brain Cortex HT, now going to identify maximums")
#
# identify_maximum_pext_per_gene(context_salmon_brain_cortex_annotated, context_salmon_brain_cortex_max_per_gene)

print("starting salmon")
mt, salmon = read_tx_annotation_tables(gnomad_release_mt_path, salmon_brain_cortex_out_path, 'ht')

mt_annotated = tx_annotate_mt(mt, salmon,"proportion", tissues_to_filter=None,gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path,filter_to_csqs=all_coding_csqs)
mt_annotated.write("gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/data/gnomad.exomes.r2.1.1.sites.tx_annotated.brain.cortex.salmon.020520.mt",overwrite = True)
