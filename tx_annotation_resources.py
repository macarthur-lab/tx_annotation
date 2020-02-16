from gnomad_hail import *

# MTs/HTs of interest
gnomad_release_mt_path = "gs://gnomad-public/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht"
gnomad_genomes_release_mt_path = "gs://gnomad-public/release/2.1.1/ht/genomes/gnomad.genomes.r2.1.1.sites.ht"
clinvar_ht_path = "gs://gnomad-resources/clinvar/hail-0.2/clinvar_20181028.vep.ht" #available in gnomad_hail as well
ddid_asd_de_novos = "gs://gnomad-public/papers/2019-tx-annotation/results/de_novo_variants/asd_ddid_de_novos.vepped.021819.mt"
chd_de_novos = "gs://gnomad-berylc/tx-annotation/hail2/DeNovoSignal/chd/congenital_heart_disease_meta.minrep.dedup.both.loftee.beta.vep.121018.mt"
context_ht_path = "gs://gnomad-public/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht"

# Precomputed HTs with pext for all possible bases (calculated from context_ht_path, script for generating : get_all_possible_snvs.py
# The first version includes all of GTEx tissues, the second is the version used in the manuscript, that has certain tissues filtered (see manuscript)
all_bases_gtex_pext = "gs://gnomad-berylc/tx-annotation/hail2/browser_integration/all.baselevel.021520.tsv.bgz"
#all_bases_gtex_pext_filtered_tissues = "gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.0218219.ht"

# GTEx files
gtex_v7_tx_summary_ht_path='gs://gnomad-public/papers/2019-tx-annotation/data/GRCH37_hg19/GTEx.V7.tx_medians.021420.ht'
gtex_v7_gene_maximums_ht_path = "gs://gnomad-public/papers/2019-tx-annotation/data/GRCH37_hg19/GTEx.v7.gene_expression_per_gene_per_tissue.021420.ht"

gtex_v8_tx_summary_mt_path = "gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.V8.tx_medians.042319.ht"
gtex_v8_gene_maximums_ht_path = "gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.v8.gene_expression_per_gene_per_tissue.042319.kt"

# Gene lists
curated_haploinsufficient_genes = "gs://gnomad-public/papers/2019-tx-annotation/data/gene_lists/HI_genes_100417.tsv"
recessive_disease_genes = 'gs://gnomad-public/papers/2019-tx-annotation/data/gene_lists/all_ar.tsv'
constraint = "gs://gnomad-resources/lof_paper/full_lof_metrics_by_transcript_an_adj_by_gene.txt.bgz"

# CSQ terms
lof_csqs = ["stop_gained","splice_donor_variant", "splice_acceptor_variant","frameshift_variant"]
missense_csqs = ["missense_variant", "inframe_insertion", "inframe_deletion"]
syn_csqs = ["synonymous_variant"]
all_coding_csqs = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT

# Misc
v7_tissues_to_drop = ["Bladder", "Brain_Spinalcord_cervicalc_1_", "Brain_Substantianigra",
                      "Cervix_Ectocervix","Cervix_Endocervix", "FallopianTube", "Kidney_Cortex",
                      "MinorSalivaryGland", "Uterus", "Ovary","Testis", "Vagina",
                      "Cells_EBV_transformedlymphocytes", "Cells_Transformedfibroblasts", "Prostate"]
phylocsf_file_path = "gs://gnomad-public/papers/2019-tx-annotation/data/other_data/phylocsf_data.tsv.bgz"


# Fetal RNA-seq
hbdr_fetal_tissue_summary_ht_path = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.tissue_names.medians.021520.mt"
hbdr_fetal_tissue_gene_maximums_ht_path ="gs://gnomad-public/papers/2019-tx-annotation/data/HBDR_fetal_RNAseq/HBDR.gene_expression_per_gene_per_tissue.021520.ht"


def make_clinvar_hail2(clinvar_vcf_path, clinvar_variants_table, clinvar_mt_out_path):
    """
    Import ClinVar vcf file, and turn it into a usable Hail2 mt

    :param str clinvar_vcf_path: Example : "gs://gnomad-berylc/tx-annotation/hail2/clinvar_alleles_single.b37.vcf.bgz"
    :param str clinvar_variants_table: Example : "gs://gnomad-berylc/tx-annotation/hail2/clinvar_alleles_single.b37.variants_table.tsv"
    :param bool repartition:
    :param int n_partitions: Number of partitions if repartition = True
    :param str clinvar_mt_out_path: "gs://gnomad-resources/clinvar/hail-0.2/clinvar_alleles.single.b37.hail2.vepped.mt"
    :return: split and VEP'd MT
    :rtype: MatrixTable
    """
    clinvar_mt = hl.import_vcf(clinvar_vcf_path)
    variants_table = hl.import_table(clinvar_variants_table, impute=True)
    variants_table = variants_table.annotate(v = hl.parse_variant(variants_table.v))
    variants_table = (variants_table.annotate(locus=variants_table.v.locus, alleles=variants_table.v.alleles).
                      key_by('locus', 'alleles'))

    clinvar_mt = clinvar_mt.annotate_rows(va=variants_table[clinvar_mt.locus, clinvar_mt.alleles])

    clinvar_mt = split_multi_dynamic(clinvar_mt, left_aligned=False)
    clinvar_mt = clinvar_mt.repartition(100)
    clinvar_vep = hl.vep(clinvar_mt, vep_config)
    clinvar_vep.write(clinvar_mt_out_path, overwrite=True)

    t = hl.read_matrix_table(clinvar_mt_out_path)
    t.rows().show()


def revep_with_loftee_beta(vepped_mt_path, out_mt_path):
    clinvar = hl.read_matrix_table(clinvar_vepped_mt_path)
    clinvar = clinvar.drop(clinvar.vep)
    clinvar_vep = hl.vep(clinvar_mt, vep_config)
    clinvar_vep.write(out_mt_path)


def vcf_to_hail2(vcf_path, mt_out):
    """
    From Konrad, who had already did this. I didn't actually run the code.

    :param str vcf_path: Example   "gs://gnomad/raw/source/ExAC.r1.sites.vep.vcf.gz"
    :param str mt_out: Example: "gs://gnomad/raw/hail-0.2/vds/exac/exac.r1.sites.vep.vds" should be mt but whatevs
    :param bool repartition:
    :param int n_partitions: Number of partitions if repartition = True
    :return: Writes out VEP'd Hail0.2 MatrixTable
    :rtype: None
    """
    mt = hl.import_vcf(vcf_path, force_bgz=True, min_partitions=1000)
    mt = split_multi_dynamic(mt)
    mt = hl.vep(mt, vep_config)
    mt.write(mt_out)

def make_gnomad_release_hail2(vcf_path, mt_out):
    """
    Used to import, filter and VEP existing "bootleg" gnomAD VCF (01.26.2018) and write out as a Hail 0.2 MatrixTable

    :param str vcf_path:
    Example: "gs://gnomad-public/release/2.0.2/vcf/exomes/gnomad.exomes.r2.0.2.sites.vcf.bgz"
    :param str mt_out:
    Example: "gs://gnomad-berylc/tx-annotation/hail2/gnomad.exomes.r2.0.2.sites.split.vep.030818.mt"
    :return: Writes out VEP'd Hail0.2 MatrixTable
    :rtype: None
    """
    release_mt = hl.import_vcf(vcf_path, min_partitions=8000)
    release_mt = split_multi_dynamic(release_mt)

    release_mt = release_mt.annotate_rows(
        as_pass=(release_mt.info.AS_FilterStatus[release_mt.a_index - 1] == "PASS") & (
                release_mt.filters.length() == 0))

    release_mt = release_mt.filter_rows(release_mt.as_pass)
    release_mt = hl.vep(release_mt, vep_config)
    release_mt.write(mt_out)
    # Confirmed 13,443,237 variants


def get_gtex_summary(gtex_rsem_path, gtex_tx_summary_out_path, get_medians=True, make_per_tissue_file=None):
    """
    Get GTEx RSEM table with ENSTs and ENSGs as rows and GTEx samples as columns (e.g. Muscle-Skeletal.12,
    Adipose.27 etc.) and write out a table with same rows, and tissues as columns (Muscle-Skeletal, Adipose etc.)
    with cells representing summary expression of transcripts across tissues (ie. mean or median).

    :param str gtex_rsem_path: Output of RSEM quantifications from GTEx
    Example: "gs://gnomad-berylc/tx-annotation/reheadered.031216.GTEx_Analysis_2016-09-07_RSEMv1.2.22_transcript_tpm.txt.bgz"
    :param str gtex_tx_summary_out_path: Path to write out.
    Example: "gs://gnomad-berylc/tx-annotation/hail2/GTEx.V7.tx_medians.030818.mt"
    :param bool get_medians: Default True. If False, returns mean transcript expression per tissue
    :return: Writes out summarized GTEx transcript expression as Table.
    :rtype: None
    """
    gtex = hl.import_matrix_table(gtex_rsem_path, row_key='transcript_id',
                                  row_fields={'transcript_id': hl.tstr, 'gene_id': hl.tstr}, entry_type=hl.tfloat64,
                                  force_bgz=True)

    gtex = gtex.annotate_cols(tissue=gtex.col_id.split("\\.")[0])

    # Will come in handy when it's time to splat
    tissue_ids = sorted(gtex.aggregate_cols(hl.agg.collect_as_set(gtex.tissue)))
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}

    if get_medians:
        gtex = gtex.group_cols_by(gtex.tissue).aggregate(tx_expr_summary=hl.median(hl.agg.collect(gtex.x)))
    else:
        gtex = gtex.group_cols_by(gtex.tissue).aggregate(tx_expr_summary=hl.mean(hl.agg.collect(gtex.x)))

    # Make a new column as an array of the values across tissues (per transcript)
    # Added sorting to ensure match up, a bit more wordy then before but oh well.
    gtex_table = gtex.entries()
    gtex_table = gtex_table.key_by(gtex_table.transcript_id, gtex_table.gene_id)
    gtex_table = gtex_table.collect_by_key()
    gtex_table = gtex_table.annotate(agg_expression=hl.sorted(gtex_table.values, key=lambda x: x.tissue).
                                     map(lambda x: x.tx_expr_summary))

    # Strip version numbers
    gtex_table = gtex_table.key_by(transcript_id=gtex_table.transcript_id.split("\\.")[0],
                                   gene_id=gtex_table.gene_id.split("\\.")[0])

    # Write out ht
    gtex_table.write(gtex_tx_summary_out_path, overwrite=True)

    if make_per_tissue_file:
        gtex_table_modified = gtex_table.annotate(**{
            tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
                gtex_table.agg_expression[d[tissue_id]] for tissue_id in tissue_ids})

        gtex_table_modified = gtex_table_modified.drop(gtex_table_modified.values)

        gtex_table_modified.export(make_per_tissue_file)



def get_gene_expression(gtex_tx_summary_path, gene_expression_out):
    gtex = hl.read_table(gtex_tx_summary_path)
    gene_expression = gtex.group_by(ensg=gtex.gene_id).aggregate(gene_expression=hl.agg.array_sum(gtex.agg_expression))
    tissue_ids = sorted([y.tissue for y in gtex.values.take(1)[0]])
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}
    gene_expression = (gene_expression.annotate(
        gene_maximum_dict={tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
                               gene_expression.gene_expression[d[tissue_id]] for tissue_id in tissue_ids}))

    gene_expression.show(10)
    gene_expression.write(gene_expression_out, overwrite=True)


def import_and_modify_gene_maximums(gtex_gene_maximums_table_tsv_path, gtex_gene_maximums_table_kt_out_path):
    """
    This function is deprecated
    :param gtex_v7_gene_maximums_table_path:
    Example: =  "gs://gnomad-berylc/tx-annotation/hail2/data/GTEx.v7.max_expression_per_base_per_tissue.021118.tsv"):
    :return:
    """
    gene_maximum_kt = hl.import_table(gtex_gene_maximums_table_tsv_path, impute=True, force_bgz=True)
    if 'Name' in list(gene_maximum_kt.row):
        gene_maximum_kt = gene_maximum_kt.rename({"Name": "ensg"})

    tissues = list(set(gene_maximum_kt.row) - set(['ensg', 'Description','representatative bases']))

    gene_maximum_kt = gene_maximum_kt.annotate(
        gene_maximum_dict={tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
                               gene_maximum_kt[tissue_id] for tissue_id in tissues})

    gene_maximum_kt = gene_maximum_kt.drop(*tissues)

    gene_maximum_kt = gene_maximum_kt.key_by(ensg=gene_maximum_kt.ensg.split("\\.")[0])


    gene_maximum_kt.write(gtex_gene_maximums_table_kt_out_path, overwrite=True)

def identify_maximum_pext_per_gene(all_possible_snv_pext_path, gene_maximums_out):
    """
    :param all_possible_snv_pext_path: Path to HT with pext annotations for all possible bases
    For manuscript, used "gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021819.ht"
    :param gene_maximums_out: Path to file that will have have the maximum pext per gene, which is used to filter genes with low maximum pexts
    For manuscript : "gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.genes.max.pext.010820.tsv.bgz"
    :return:
    """
    ht = hl.read_matrix_table(all_possible_snv_pext_path)
    ht = ht.filter_rows(~hl.is_missing(ht.tx_annotation))

    # Filter only to variants that have a tx annotation value (aka all coding/exonic variants)
    ht = ht.filter_rows(~hl.is_missing(ht.tx_annotation))

    # Spread out all values for a variant
    ht_exploded = ht.explode_rows(ht.tx_annotation)
    ht_exploded = ht_exploded.annotate_rows(
        all_tissues=ht_exploded.tx_annotation.drop('ensg', 'csq', 'symbol', 'lof', 'lof_flag'))

    # For each gene, get the maximum pext value - we are interested in genes where the max pext value is low (0.2 cutoff for manuscript)
    max_per_gene = ht_exploded.group_rows_by(
        ht_exploded.tx_annotation.ensg, ht_exploded.tx_annotation.symbol).aggregate_rows(
        max_pexts=hl.struct(**{key: hl.agg.max(val) for key, val in ht_exploded.all_tissues.items()})).result()

    max_per_gene = max_per_gene.annotate_rows(**max_per_gene.max_pexts)

    max_per_gene.rows().export(gene_maximums_out)

