from tx_annotation import *
mt = get_gnomad_data('exomes', adj=True, non_refs_only=True,release_samples=True, release_annotations=True)
hi_variants = hl.import_table("gs://gnomad-berylc/tx-annotation/gnomad_release/ellie_curation/gnomad.HIgene.variants.for.curation.113018.tsv")

# Filter variants

rt = hl.import_table("gs://gnomad-berylc/tx-annotation/gnomad_release/ellie_curation/gnomad.HIgene.variants.for.curation.113018.tsv")
rt = rt.annotate(variant=rt.chrom + ':' + rt.pos + ":" + rt.ref + ":" + rt.alt)
rt = rt.annotate(annotation = "ellie")
rt = rt.annotate(** hl.parse_variant(rt.variant))

mt =  mt.annotate_rows(in_curation = rt[mt.locus, mt.alleles])


# Annotate
vep_ht = hl.read_table(annotations_ht_path('exomes', 'vep_csq'))

mt = mt.annotate_rows(info=hl.struct(
    AC=mt.freq[0].AC, AN=mt.freq[0].AN, AF=mt.freq[0].AF, n_hom=mt.freq[0].homozygote_count,
    CSQ=vep_ht[mt.row_key].vep)).drop('is_missing')

hl.export_vcf(mt,"gs://gnomad-berylc/tx-annotation/gnomad_release/ellie_curation/gnomad.HIgene.variants.for.curation.011418.vcf.bgz",
              metadata={'info': {'CSQ': {'Description': vep_ht.vep_csq_header.collect()[0]}}})