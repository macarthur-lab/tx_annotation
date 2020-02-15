'''Code snippet to annotate all possible SNVs pext values. Ran on 100 non-preemptible nodes'''

from tx_annotation import *
from datetime import datetime
'''

print(gtex_v7_tx_summary_ht_path)
print(gtex_v7_gene_maximums_ht_path)

start = datetime.now()
print("Starting on",start)

mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_ht_path, 'ht')
start_annotation = datetime.now()

print("Starting tx annotation on", start_annotation)

context_hg19_annotated = tx_annotate_mt(mt, gtex, "proportion",
                                        tissues_to_filter=None,
                                        gene_maximums_ht_path=gtex_v7_gene_maximums_ht_path,
                                        filter_to_csqs=all_coding_csqs)

finished_annotation = datetime.now()
print("Finished annotation ",finished_annotation)


context_hg19_annotated.write("gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht",overwrite=True)

# get maximum pext per gene to identify genes to filter from future analyses
identify_maximum_pext_per_gene("gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021520.ht",
                               "gs://gnomad-public/papers/2019-tx-annotation/data/GRCH37_hg19/all.genes.max.pext.GTEx.v7.021520.tsv.bgz")

stop = datetime.now()
print("Finished writing",stop)
print("total time", stop - start)

# Now post-processing to enable users to extract tx annotation without hail
start_post = datetime.now()

print("Starting post-processing for those who don't want to use Hail",start_post)

context_hg19_annotated = context_hg19_annotated.filter_rows(~hl.is_missing(context_hg19_annotated.tx_annotation))
context_hg19_annotated = context_hg19_annotated.key_rows_by()
context_hg19_annotated = context_hg19_annotated.select_rows("locus", "alleles", "tx_annotation")
context_hg19_annotated = context_hg19_annotated.annotate_rows(chrom=context_hg19_annotated.locus.contig,
                                                              pos=context_hg19_annotated.locus.position,
                                                              ref=context_hg19_annotated.alleles[0],
                                                              alt=context_hg19_annotated.alleles[1])

context_hg19_annotated = context_hg19_annotated.select_rows("chrom", "pos", "ref", "alt", "tx_annotation")
context_hg19_annotated.count_rows()
context_hg19_annotated.rows().export("gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.GTEx.v7.021520.tsv.bgz")
'''

