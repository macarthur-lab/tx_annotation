'''
Analysis to get a sense of whether pext is robust to isoform quantification tool used

Because it would be difficult to reprocess all of GTEx, we're only using GTEx Brain Cortex samples here for comparison

This is also a good way to walk through how to prepare your own isoform expression matrix for pext/ext annotation

I have two isoform expression matrices on the same 151 Brain cortex samples. First two fields are transcript_id and
gene_id
'''
'''
# Make median tissue expression files to be used for analysis

salmon_brain_cortex_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/salmon.GTEx.v7.brain.cortex.021419.tsv"
salmon_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/salmon.GTEx.v7.brain.cortex.tx_medians.021419.mt"
get_gtex_summary(salmon_brain_cortex_path, salmon_brain_cortex_out_path, get_medians=True, make_per_tissue_file=None)

rsem_brain_cortex_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.GTEx.v7.brain.cortex.021419.tsv"
rsem_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.GTEx.v7.brain.cortex.tx_medians.021419.mt"
get_gtex_summary(rsem_brain_cortex_path, rsem_brain_cortex_out_path, get_medians=True, make_per_tissue_file=None)

# Get baselevel TCF4 pext values on both
salmon_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/salmon.GTEx.v7.brain.cortex.tx_medians.021519.mt"
mt, salmon = read_tx_annotation_tables(context_ht_path, salmon_brain_cortex_out_path, "ht")
salmon_baselevel_prop= get_baselevel_expression_for_genes(mt, salmon, gene_list = {'TCF4'}, get_proportions=True)
salmon_baselevel_prop.export("gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/TCF4/salmon.TCF4.baselevel.pext.021519.tsv.bgz")

rsem_brain_cortex_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.GTEx.v7.brain.cortex.tx_medians.021519.mt"
mt, rsem = read_tx_annotation_tables(context_ht_path, rsem_brain_cortex_path, "ht")
rsem_baselevel_prop= get_baselevel_expression_for_genes(mt, rsem, gene_list = {'TCF4'}, get_proportions=True)
rsem_baselevel_prop.export("gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/TCF4/rsem.TCF4.baselevel.pext.021519.tsv.bgz")

#  To run MAPs on both, first annotate all of gnomAD exomes with both files, then follow the maps_submit_per_class.py script
# replacing the HT on line 32 with the files you create here
# First had to annotate gnomAD with each e.g.: 
latest_gnomad_120518 = "gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht"
rsem_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.GTEx.v7.brain.cortex.tx_medians.021519.mt"
mt, gtex = read_tx_annotation_tables(latest_gnomad_120518, rsem_brain_cortex_out_path, 'ht')
mt_annotated = tx_annotate_mt(mt, gtex,"proportion", tissues_to_filter = None,
                              gene_maximums_ht_path = gtex_v7_gene_maximums_ht_path,
                              filter_to_csqs=all_coding_csqs)
mt_annotated.write("gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.brain_cortex.gnomad.exomes.r2.1.1.sites.tx_annotated.021519.ht")
'''

# Gene list annotation

# Salmon
# gnomAD
salmon_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/salmon.GTEx.v7.brain.cortex.tx_medians.021519.mt"
clinvar_mt, salmon = read_tx_annotation_tables(clinvar_ht_path, salmon_brain_cortex_out_path, "ht")
mt, salmon = read_tx_annotation_tables(gnomad_release_mt_path, salmon_brain_cortex_out_path, "ht")
mt = mt.filter_rows(hl.len(mt.filters) == 0)

hi_genes = import_gene_list(curated_haploinsufficient_genes, gene_column="ENSGID", ensg=True)

out_dir = "gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/gene_list_comparisons_salmon_rsem/"


mt_gnomad_hi = tx_annotate_mt(mt, salmon,"proportion", tissues_to_filter = None,
                              filter_to_csqs=lof_csqs,
                              filter_to_genes=hi_genes, gene_column_in_mt="gene_id")
mt_gnomad_hi = mt_gnomad_hi.filter_rows(~hl.is_missing(mt_gnomad_hi.tx_annotation))
mt_gnomad_hi = pull_out_worst_from_tx_annotate(mt_gnomad_hi)
mt_gnomad_hi.rows().export("%sHI_genes.salmon.brain_cortex.gnomad.exomes.r2.1.tx_annotated.021520.tsv.bgz" %out_dir)

# ClinVar
mt_clinvar_hi = tx_annotate_mt(clinvar_mt, salmon,"proportion",tissues_to_filter = None,
                               filter_to_csqs=lof_csqs, filter_to_genes=hi_genes,
                               gene_column_in_mt="gene_id")
mt_clinvar_hi = mt_clinvar_hi.filter_rows(~hl.is_missing(mt_clinvar_hi.tx_annotation))
mt_clinvar_hi = pull_out_worst_from_tx_annotate(mt_clinvar_hi)
mt_clinvar_hi = mt_clinvar_hi.annotate_rows(**mt_clinvar_hi.info)
mt_clinvar_hi = mt_clinvar_hi.drop("vep", "tx_annotation","info")
mt_clinvar_hi.rows().export("%sHI_genes.salmon.brain_cortex.clinvar.alleles.single.b37.tx_annotated.021520.tsv.bgz" %out_dir)


# RSEM
# gnomAD
rsem_brain_cortex_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/rsem.GTEx.v7.brain.cortex.tx_medians.021519.mt"
clinvar_mt, rsem = read_tx_annotation_tables(clinvar_ht_path, rsem_brain_cortex_path, "ht")
mt, rsem = read_tx_annotation_tables(gnomad_release_mt_path, rsem_brain_cortex_path, "ht")
mt = mt.filter_rows(hl.len(mt.filters) == 0)

hi_genes = import_gene_list(curated_haploinsufficient_genes, gene_column="ENSGID", ensg=True)

out_dir = "gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/gene_list_comparisons_salmon_rsem/"


mt_gnomad_hi = tx_annotate_mt(mt, rsem,"proportion", tissues_to_filter = None,
                              filter_to_csqs=lof_csqs,
                              filter_to_genes=hi_genes, gene_column_in_mt="gene_id")
mt_gnomad_hi = mt_gnomad_hi.filter_rows(~hl.is_missing(mt_gnomad_hi.tx_annotation))
mt_gnomad_hi = pull_out_worst_from_tx_annotate(mt_gnomad_hi)
mt_gnomad_hi.rows().export("%sHI_genes.rsem.brain_cortex.gnomad.exomes.r2.1.tx_annotated.021520.tsv.bgz" %out_dir)

mt_clinvar_hi = tx_annotate_mt(clinvar_mt, rsem,"proportion",tissues_to_filter = None,
                               filter_to_csqs=lof_csqs, filter_to_genes=hi_genes,
                               gene_column_in_mt="gene_id")
mt_clinvar_hi = mt_clinvar_hi.filter_rows(~hl.is_missing(mt_clinvar_hi.tx_annotation))
mt_clinvar_hi = pull_out_worst_from_tx_annotate(mt_clinvar_hi)
mt_clinvar_hi = mt_clinvar_hi.annotate_rows(**mt_clinvar_hi.info)
mt_clinvar_hi = mt_clinvar_hi.drop("vep", "tx_annotation","info")
mt_clinvar_hi.rows().export("%sHI_genes.rsem.brain_cortex.clinvar.alleles.single.b37.tx_annotated.021520.tsv.bgz" %out_dir)

# Don't forget gnomAD genomes!
'''
salmon_brain_cortex_out_path = "gs://gnomad-berylc/tx-annotation/hail2/salmon_rsem/salmon.GTEx.v7.brain.cortex.tx_medians.021520.mt"
clinvar_mt, salmon = read_tx_annotation_tables(clinvar_ht_path, salmon_brain_cortex_out_path, "ht")
mt, salmon = read_tx_annotation_tables(gnomad_genomes_release_mt_path, salmon_brain_cortex_out_path, "ht")
mt = mt.filter_rows(hl.len(mt.filters) == 0)

hi_genes = import_gene_list(curated_haploinsufficient_genes, gene_column="ENSGID", ensg=True)

out_dir = "gs://gnomad-public/papers/2019-tx-annotation/results/salmon_rsem/gene_list_comparisons_salmon_rsem/"


mt_gnomad_hi = tx_annotate_mt(mt, salmon,"proportion", tissues_to_filter =None,
                              filter_to_csqs=lof_csqs,
                              filter_to_genes=hi_genes, gene_column_in_mt="gene_id")
mt_gnomad_hi = mt_gnomad_hi.filter_rows(~hl.is_missing(mt_gnomad_hi.tx_annotation))
mt_gnomad_hi = pull_out_worst_from_tx_annotate(mt_gnomad_hi)
mt_gnomad_hi.rows().export("%sHI_genes.salmon.brain_cortex.gnomad.genomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)
'''