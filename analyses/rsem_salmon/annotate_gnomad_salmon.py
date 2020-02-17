'''Code snippet to annotate gnomAD sites MT with pext values. Ran on 50 non-preemptible nodes'''

from tx_annotation import *
from datetime import datetime

start = datetime.now()
print("Starting on",start)
latest_gnomad_120518 = "gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht"

base = "gs://gnomad-berylc/tx-annotation/hail2/reviewer_response/maximum_0/salmon_requantification_pc_only/"
# salmon_brain_pc_only = "%srenamed.updated.salmon.GTEx.v7.brain.cortex.pc_transcripts.only.medians.ht"%base
# salmon_brain_pc_only_max = "%srenamed.updated.salmon.GTEx.v7.brain.cortex.pc_transcripts.only.maximums.ht"%base
#mt, salmon_pc = read_tx_annotation_tables(latest_gnomad_120518, salmon_brain_pc_only, 'ht')

base = "gs://gnomad-berylc/tx-annotation/hail2/reviewer_response/maximum_0/salmon_requantification_pc_only/"
salmon_brain_lncrna = "%srenamed.salmon.GTEx.v7.brain.cortex.pc_lncRNA.medians.ht"%base
salmon_brain_lncrna_max = "%srenamed.salmon.GTEx.v7.brain.cortex.pc_lncRNA.maximums.ht"%base
mt, salmon_lncrna = read_tx_annotation_tables(latest_gnomad_120518, salmon_brain_lncrna, 'ht')


start_annotation = datetime.now()
print("Starting tx annotation on",start_annotation)
# mt_annotated = tx_annotate_mt(mt, salmon_pc,"proportion",
#                               gene_maximums_kt_path=salmon_brain_pc_only_max,
#                               filter_to_csqs=all_coding_csqs,
#                               tissues_to_filter=None)

mt_annotated = tx_annotate_mt(mt, salmon_lncrna,"proportion",
                              gene_maximums_kt_path=salmon_brain_lncrna_max,
                              filter_to_csqs=all_coding_csqs,
                              tissues_to_filter=None)

finished_annotation = datetime.now()
print("Finished annotation ",finished_annotation)

# mt_annotated.write("%sgnomad.exomes.r2.1.1.sites.salmon_pc_only.083119.ht"%base,overwrite = True)
mt_annotated.write("%sgnomad.exomes.r2.1.1.sites.salmon_pc_lncRNA.090119.ht"%base,overwrite = True)

print("Done")

