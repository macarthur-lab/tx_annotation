'''Code snippet to annotate all possible SNVs pext values. Ran on 100 non-preemptible nodes'''

from tx_annotation import *
from datetime import datetime

baselevel = hl.read_matrix_table("gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.022719.ht")
baselevel = baselevel.filter_rows(~hl.is_missing(baselevel.tx_annotation))
baselevel = baselevel.key_rows_by()
baselevel = baselevel.select_rows("locus", "alleles", "tx_annotation")
baselevel = baselevel.annotate_rows(chrom = baselevel.locus.contig,
                                    pos = baselevel.locus.position,
                                    ref = baselevel.alleles[0],
                                    alt = baselevel.alleles[1])
baselevel = baselevel.select_rows("chrom", "pos", "ref", "alt", "tx_annotation")
baselevel.count_rows()
baselevel.rows().export("gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.022719.tsv.bgz")

'''

start = datetime.now()
print("Starting on",start)

mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_mt_path, 'ht')

start_annotation = datetime.now()
print("Starting tx annotation on", start_annotation)
mt_annotated = tx_annotate_mt(mt, gtex, "proportion",
                              tissues_to_filter=None,
                              gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path,
                              filter_to_csqs=all_coding_csqs)
finished_annotation = datetime.now()
print("Finished annotation ",finished_annotation)

mt_annotated.write("gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.022719.ht", overwrite=True)

#mt_annotated.rows().export("gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021819.tsv.bgz")

stop = datetime.now()
print("Finished writing",stop)
print("total time", stop - start)

# Now post-processing and exporting to TSV
start_post = datetime.now()
print("Starting post-processing for those who don't want to use Hail",start_post)

mt_annotated = mt_annotated.filter_rows(~hl.is_missing(mt_annotated.tx_annotation))
mt_annotated = mt_annotated.annotate_rows(tx_annotation=mt_annotated.tx_annotation.map(fix_loftee_beta_nonlofs))
mt_annotated = pull_out_worst_from_tx_annotate(mt_annotated)

mt_annotated = mt_annotated.select_rows('ensg' ,
                                        'csq',
                                        'symbol',
                                        'lof'  ,
                                        'lof_flag',
                                        'Cells_Transformedfibroblasts',
                                        'Prostate',
                                        'Spleen'  ,
                                        'Brain_FrontalCortex_BA9_'  ,
                                        'SmallIntestine_TerminalIleum'  ,
                                        'MinorSalivaryGland',
                                        'Artery_Coronary'  ,
                                        'Skin_SunExposed_Lowerleg_',
                                        'Cells_EBV_transformedlymphocytes',
                                        'Brain_Hippocampus' ,
                                        'Esophagus_Muscularis',
                                        'Brain_Nucleusaccumbens_basalganglia_'  ,
                                        'Artery_Tibial'  ,
                                        'Brain_Hypothalamus'  ,
                                        'Adipose_Visceral_Omentum_'  ,
                                        'Cervix_Ectocervix',
                                        'Brain_Spinalcord_cervicalc_1_',
                                        'Brain_CerebellarHemisphere'  ,
                                        'Nerve_Tibial',
                                        'Breast_MammaryTissue',
                                        'Liver'  ,
                                        'Skin_NotSunExposed_Suprapubic_'  ,
                                        'AdrenalGland'  ,
                                        'Vagina',
                                        'Pancreas'  ,
                                        'Lung'  ,
                                        'FallopianTube',
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
                                        'Brain_Substantianigra',
                                        'Kidney_Cortex',
                                        'Uterus',
                                        'Stomach'  ,
                                        'WholeBlood',
                                        'Bladder',
                                        'Brain_Anteriorcingulatecortex_BA24_',
                                        'Brain_Putamen_basalganglia_'  ,
                                        'Brain_Caudate_basalganglia_'  ,
                                        'Colon_Sigmoid'  ,
                                        'Cervix_Endocervix',
                                        'Ovary',
                                        'Esophagus_GastroesophagealJunction'  ,
                                        'Brain_Amygdala',
                                        'Testis',
                                        'mean_proportion')

mt_annotated.rows().export("gs://gnomad-public/papers/2019-tx-annotation/data/gnomad.exomes.r2.1.1.sites.tx_annotated.021319.tsv.bgz")

stop_post = datetime.now()
print("Finished post-processing for figures", stop_post)


'''