from tx_annotation import *
'''
rt = hl.import_table("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/updated_includeASC_Hamdam/ASD_ddid_Updated.hamdam.ASC.011418.txt")
rt = rt.annotate(variant=rt.CHROM + ':' + rt.POSITION + ":" + rt.REF + ":" + rt.ALT)
rt = rt.annotate(** hl.parse_variant(rt.variant))
rt = rt.key_by(rt.locus, rt.alleles)
mt = hl.MatrixTable.from_rows_table(rt)
mt = mt.repartition(10)
annotated_mt = hl.vep(mt, vep_config)
annotated_mt.write("gs://gnomad-public/papers/2019-tx-annotation/results/de_novo_variants/asd_ddid_de_novos.vepped.021819.mt")
'''


mt, gtex = read_tx_annotation_tables(ddid_asd_de_novos, gtex_v7_tx_summary_mt_path, "mt")
ddid_asd = tx_annotate_mt(mt, gtex, "proportion",
                          filter_to_csqs=all_coding_csqs)
ddid_asd = ddid_asd.filter_rows(~hl.is_missing(ddid_asd.tx_annotation))
ddid_asd = ddid_asd.annotate_rows(tx_annotation=ddid_asd.tx_annotation.map(fix_loftee_beta_nonlofs))
ddid_asd = pull_out_worst_from_tx_annotate(ddid_asd)
ddid_asd.rows().export("gs://gnomad-public/papers/2019-tx-annotation/results/de_novo_variants/"
                       "asd_ddid_de_novos.tx_annotated.021819.tsv.bgz")