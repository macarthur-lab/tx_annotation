

from tx_annotation import *
rt = hl.import_table("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/all.denovo.variants.032518.txt")
rt = rt.annotate(variant=rt.CHROM + ':' + rt.POSITION + ":" + rt.REF + ":" + rt.ALT)
rt = rt.annotate(** hl.parse_variant(rt.variant))
rt = rt.key_by(rt.locus, rt.alleles)
mt = hl.MatrixTable.from_rows_table(rt)
mt = mt.repartition(10)
annotated_mt = hl.vep(mt, vep_config)
annotated_mt.write("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/all.denovo.variants.vepped.032518.mt")
#LOFTEE BETA
annotated_mt.write("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/all.denovo.variants.loftee.beta.vepped.100818.mt")


mt, gtex = read_tx_annotation_tables("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/all.denovo.variants.vepped.032518.mt",  gtex_v7_tx_summary_mt_path)
tx_table = tx_annotate_mt(mt, gtex, "proportion", filter_to_csqs=all_coding_csqs)



# with new annotation
mt, gtex = read_tx_annotation_tables("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/all.denovo.variants.vepped.032518.mt",  gtex_v7_tx_summary_mt_path)

mt_annotated = mt_annotated.annotate_rows(**hl.sorted( mt_annotated.tx_annotation, key = lambda x : csq_order[(x.lof, x.lof_flag == "",x.csq)])[0])


mt_annotated.rows().export("gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/data/merged.in.hail.all.denovo.variants.vepped.annotated.proportions.100118.tsv.bgz")
