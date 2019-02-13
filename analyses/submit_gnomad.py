'''Code snippet to annotate gnomAD sites MT with pext values. Ran on 50 non-preemptible nodes'''

from tx_annotation import *
from datetime import datetime

start = datetime.now()
print("Starting on",start)
latest_gnomad_120518 = "gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht"

mt, gtex = read_tx_annotation_tables(latest_gnomad_120518, gtex_v7_tx_summary_mt_path, 'ht')

start_annotation = datetime.now()
print("Starting tx annotation on",start_annotation)
mt_annotated = tx_annotate_mt(mt, gtex,"proportion", gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path,filter_to_csqs=all_coding_csqs)
finished_annotation = datetime.now()
print("Finished annotation ",finished_annotation)

mt_annotated.write("gs://gnomad-public/papers/2019-tx-annotation/gnomad.exomes.r2.1.1.sites.tx_annotated.021319.ht",overwrite = True)

stop = datetime.now()
print("Finished writing",stop)
print("total time", stop - start)

# Now post-processing and exporting to TSV
start_post = datetime.now()
print("Starting post-processing for figures",start_post)

mt_annotated = mt_annotated.filter_rows(~hl.is_missing(mt_annotated.tx_annotation))
mt_annotated = mt_annotated.annotate_rows(tx_annotation=mt_annotated.tx_annotation.map(fix_loftee_beta_nonlofs))
mt_annotated = pull_out_worst_from_tx_annotate(mt_annotated)

mt_annotated = mt_annotated.select_rows('ensg' ,
                                        'csq',
                                        'symbol',
                                        'lof'  ,
                                        'lof_flag',
                                        'Spleen'  ,
                                        'Brain_FrontalCortex_BA9_'  ,
                                        'SmallIntestine_TerminalIleum'  ,
                                        'Skin_SunExposed_Lowerleg_'  ,
                                        'Artery_Coronary'  ,
                                        'Brain_Hippocampus' ,
                                        'Esophagus_Muscularis',
                                        'Brain_Nucleusaccumbens_basalganglia_'  ,
                                        'Artery_Tibial'  ,
                                        'Brain_Hypothalamus'  ,
                                        'Adipose_Visceral_Omentum_'  ,
                                        'Nerve_Tibial'  ,
                                        'Brain_CerebellarHemisphere'  ,
                                        'Liver'  ,
                                        'Breast_MammaryTissue'  ,
                                        'Skin_NotSunExposed_Suprapubic_'  ,
                                        'AdrenalGland'  ,
                                        'Pancreas'  ,
                                        'Lung'  ,
                                        'Pituitary'  ,
                                        'Muscle_Skeletal'  ,
                                        'Colon_Transverse'  ,
                                        'Artery_Aorta'  ,
                                        'Heart_AtrialAppendage'  ,
                                        'Adipose_Subcutaneous'  ,
                                        'Esophagus_Mucosa'  ,
                                        'Heart_LeftVentricle',
                                        'Brain_Cerebellum'  ,
                                        'Brain_Cortex'  ,
                                        'Thyroid'  ,
                                        'Stomach'  ,
                                        'WholeBlood',
                                        'Brain_Putamen_basalganglia_'  ,
                                        'Brain_Anteriorcingulatecortex_BA24_'  ,
                                        'Brain_Caudate_basalganglia_'  ,
                                        'Colon_Sigmoid'  ,
                                        'Esophagus_GastroesophagealJunction'  ,
                                        'Brain_Amygdala',
                                        'mean_proportion')

mt_annotated.rows().export("gs://gnomad-public/papers/2019-tx-annotation/data/gnomad.exomes.r2.1.1.sites.tx_annotated.021319.tsv.bgz")

stop_post = datetime.now()
print("Finished post-processing for figures", stop_post)
