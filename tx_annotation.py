import hail as hl
from .tx_annotation_resources import *
hl.init()


def import_gene_list(gene_list_path, gene_column, ensg=False, oe_threshold=False, peek=False):
    """
    Imports a gene list tsv and returns a set of ENSG or gene symbols

    :param str gene_list_path: Path to TSV file with gene list of interest
    :param str or None gene_column: Column in TSV file that specifies gene symbol or ENSG id.
    This column will be turned into a set.
    :param str or bool ensg: If there are no ENSGs with version numbers in the file, specify False (Default)
    If there are ENSGs with version numbers in the file, specify column containing the ENSGs.
    :param float or bool pLI_threshold: If the file does not contain pLI scores to filter, specify False (Default)
    If the file contains pLI scores, specify threshold to filter files.
    e.g. pLI threshold = 0.95
    :param bool peek: Default False.
    If you want to peek at the gene list to get the parameters
    Print out the first few lines of the gene list tsv, returns None
    :return: Set of genes of interest
    :rtype: set or None
    """
    genes = hl.import_table(gene_list_path, impute=True)
    if peek:
        genes.show(width=200)
        return None

    if oe_threshold:
        genes = genes.filter(genes.oe_lof_upper < oe_threshold)

    if ensg:
        genes = genes.annotate(ensg=genes[gene_column].split("\\.")[0])
        gene_column = "ensg"

    genes = genes.aggregate(hl.agg.collect_as_set(genes[gene_column]))

    return genes


def filter_table_to_gene_list(mt_kt, genes, gene_column_in_mt_kt):
    """Take a matrix table and return a table filtered down to a set of genes

    :param Table mt_kt:
    :param list of str or set of str genes: Genes of interest to which to filter table
    :param str gene_column_in_mt_kt: Column in matrix table that contains gene information within
    vep.transcript_consequences. often ["gene_id", "gene_symbol"]
    :return: Filtered table
    :rtype: Table
    """
    gene_names = hl.literal(genes)

    mt_kt = mt_kt.annotate(
        in_gene_of_interest=gene_names.find(lambda x: mt_kt.vep.transcript_consequences[gene_column_in_mt_kt] == x))

    mt_kt = mt_kt.filter(mt_kt.in_gene_of_interest != "NA")

    return mt_kt


def filter_table_to_csqs(mt_kt, csqs):
    """Take a matrix table and return a table filtered down to a set of CSQs
    :param Table mt_kt:
    :param list of str or set of str csqs: CSQs of interest to which to filter table
    :return: Filtered matrix table
    :rtype: Table
    """
    csqs = hl.literal(csqs)

    mt_kt = mt_kt.annotate(
        in_csq_of_interest=csqs.find(lambda x: mt_kt.vep.transcript_consequences.most_severe_consequence == x))
    mt_kt = mt_kt.filter(mt_kt.in_csq_of_interest != "NA")
    return mt_kt


def filter_clinvar_to_gene_list(mt_kt, genes, gene_column_in_mt_kt):

    gene_names = hl.literal(genes)

    mt_kt = mt_kt.annotate(
        in_gene_of_interest=gene_names.find(lambda x: mt_kt[gene_column_in_mt_kt] == x))

    mt_kt = mt_kt.filter(mt_kt.in_gene_of_interest != "NA")

    return mt_kt


def read_tx_annotation_tables(mt_path, gtex_tx_summary_path, mt_type="mt"):
    if mt_type == "mt":
        mt = hl.read_matrix_table(mt_path)

    elif mt_type == "ht":
        mt = hl.read_table(mt_path)
        mt = hl.MatrixTable.from_rows_table(mt)

    gtex = hl.read_table(gtex_tx_summary_path)
    return mt, gtex


def tx_annotate_mt(mt, gtex, tx_annotation_type,
                   tissues_to_filter = v7_tissues_to_drop, gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path,
                   filter_to_csqs=all_coding_csqs, filter_to_genes=None, gene_column_in_mt=None, filter_to_homs=False,
                   out_tx_annotation_tsv=None, out_tx_annotation_kt=None):


    """
    Annotate variants in the input MatrixTable with transcript-based expression values accross GTEx. Returns Table.

    :param MatrixTable mt: Input variant file
    :param MatrixTable gtex: Input GTEx summary MatrixTable, must have transcript_id column to key by
    :param str tx_annotation_type: One of ["expression", "proportion"]. Select proportion if you'd like the
    tx_annotation values to be normalized by max expression of the gene
    :param None or list filter_to_csqs: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of consequence terms.
    :param str gene_column_in_mt: Must be set if filter_to_genes != None.
    Column in matrix table that contains gene information within vep.transcript_consequences.
    often ["gene_id", "gene_symbol"]
    :param None or list filter_to_csqs: Default None. If you'd like to filter the mt before annotating
    (decreases time) feed in a list or set of consequence terms.
    Example = ["stop_gained","splice_donor_variant", "splice_acceptor_variant","frameshift_variant"]
    :param None or str out_tx_annotation_tsv: Default None.
    If you'd like to write out the results table as a tsv, provide a tsv path
    :param None or str out_tx_annotation_kt: Default None.
    If you'd like to write out the results table as a Hail 0.2 table, provide a .kt path
    :param bool filter_to_homs: Default False
    If True, filter to variants with at least one homozygote in dataset
    :return: Table with columns: variant, worst_csq, ensg, LOFTEE LOF, LOFTEE LOF Flag, transcript-aware expression
    by GTEx Tissue
    :rtype: Table with variants annotated with transcript-aware tissue expression
    """

    #check_inputs(**locals())

    gtex_table = gtex.key_by("transcript_id")

    #mt = process_consequences(mt, penalize_flags=False)
    mt_exploded = mt.distinct_by_row()
    mt_exploded = mt_exploded.annotate_rows(vep=mt_exploded.vep.annotate(
        transcript_consequences=mt_exploded.vep.transcript_consequences.map(add_most_severe_consequence_to_consequence)))

    # Explode the mt for the transcript consequences to be able to key by transcript ID
    mt_exploded = mt_exploded.explode_rows(mt_exploded.vep.transcript_consequences)

    mt_kt = mt_exploded.rows()
    # Currently testing removal of protein coding transcripts
    mt_kt = mt_kt.filter(mt_kt.vep.transcript_consequences.biotype == "protein_coding")

    if filter_to_genes:
        print("Filtering to genes of interest")
        mt_kt = filter_table_to_gene_list(mt_kt, filter_to_genes, gene_column_in_mt)

    if filter_to_csqs:
        print("Filtering to csqs in %s" % (",".join(filter_to_csqs)))
        mt_kt = filter_table_to_csqs(mt_kt, filter_to_csqs)

    if filter_to_homs:
        print("Filtering to variants with at least 1 homozygote sample in dataset")
        #mt_kt = mt_kt.filter(mt_kt.info.Hom[mt_kt.a_index - 1] > 0)
        idx = mt_kt.globals.freq_index_dict['gnomad']
        mt_kt = mt_kt.filter(mt_kt.freq[idx].homozygote_count >= 1)

    # Annotate mt with the gtex values (ie. join them)
    mt_kt = mt_kt.annotate(tx_data=gtex_table[mt_kt.vep.transcript_consequences.transcript_id])

    # Group by gene, worst_csq and variant, and do a pairwise-sum
    grouped_table = (
        mt_kt.group_by(csq=mt_kt.vep.transcript_consequences.most_severe_consequence,
                       ensg=mt_kt.vep.transcript_consequences.gene_id,
                       symbol=mt_kt.vep.transcript_consequences.gene_symbol,
                       locus=mt_kt.locus,
                       alleles=mt_kt.alleles,
                       lof=mt_kt.vep.transcript_consequences.lof,
                       lof_flag=mt_kt.vep.transcript_consequences.lof_flags).aggregate(tx_annotation=hl.agg.array_sum(mt_kt.tx_data.agg_expression)))

    # Expand the columns from the arrays and add tissues as headers
    #tissue_ids = gtex.tissue.collect()
    # Since gtex no longer has .tissue just a new way to do this, i probably want to save it as a global at some point
    tissue_ids = sorted([y.tissue for y in gtex.values.take(1)[0]])
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}

    tx_annotation_table = grouped_table.annotate(
        **{tissue_id.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_"):
               grouped_table.tx_annotation[d[tissue_id]] for tissue_id in tissue_ids})

    tx_annotation_table = tx_annotation_table.drop(tx_annotation_table.tx_annotation)

    # First of all do you want proportions or expression?
    if tx_annotation_type == "proportion":
        print("Returning expression proportion")
        gene_maximums_kt = hl.read_table(gene_maximums_kt_path)
        tx_annotation_table = get_expression_proportion(tx_annotation_table, tissues_to_filter, gene_maximums_kt)

    #You can write the output that is exploded by variants-ensg-csq-symbol-LOFTEE-LOFTEEflag
    # and has a value for each tissue as column, either as a TSV or a KT

    if out_tx_annotation_tsv:
        print("Writing tsv file to %s" %out_tx_annotation_tsv)
        tx_annotation_table.export(out_tx_annotation_tsv)

    if out_tx_annotation_kt:
        print("Writing Table to %s" % out_tx_annotation_kt)
        tx_annotation_table.write(out_tx_annotation_kt)

    tx_annotation_table = tx_annotation_table.key_by(tx_annotation_table.locus, tx_annotation_table.alleles)
    tx_annotation_table = tx_annotation_table.collect_by_key('tx_annotation')
    mt = mt.annotate_rows(**tx_annotation_table[mt.locus, mt.alleles])

    return mt


def get_expression_proportion(tx_table, tissues_to_filter, gene_maximum_kt):

    if tissues_to_filter:
        print("Filtering tissues:", tissues_to_filter)
        tx_table = tx_table.drop(*tissues_to_filter)

    remaining_tissue_columns = list(
        set(tx_table.row) - {'locus', 'alleles','csq', 'ensg', 'symbol','lof', 'lof_flag'})

    tx_table = tx_table.annotate(
        tx_expression=
        {tissue_id: tx_table[tissue_id] for tissue_id in remaining_tissue_columns})

    tx_table = tx_table.key_by('ensg').join(gene_maximum_kt.key_by("ensg"))

    expression_proportion_table = tx_table.annotate(
        expression_proportion_dict=
        {tissue_id: tx_table.tx_expression[tissue_id] / tx_table.gene_maximum_dict[tissue_id] for tissue_id in remaining_tissue_columns})

    columns_to_drop = list(set(expression_proportion_table.row) - {'locus', 'alleles','csq', 'ensg',
                                                                   'symbol','lof', 'lof_flag','expression_proportion_dict'})

    expression_proportion_table = expression_proportion_table.drop(*columns_to_drop)

    expression_proportion_table = expression_proportion_table.annotate(
        **{tissue_id: expression_proportion_table.expression_proportion_dict[tissue_id] for tissue_id in remaining_tissue_columns})

    expression_proportion_table = expression_proportion_table.annotate(
        mean_proportion=hl.mean( hl.filter(
            lambda e : ~hl.is_nan(e),  [expression_proportion_table[tissue_id] for tissue_id in remaining_tissue_columns]) filter_missing=True))

    expression_proportion_table = expression_proportion_table.drop(
        expression_proportion_table.expression_proportion_dict).key_by(
        'locus', 'alleles', 'ensg')

    return expression_proportion_table


def pull_out_worst_from_tx_annotate(mt):
    csq_order = []
    for loftee_filter in ["HC", "LC"]:
        for no_flag in [True, False]:
            for consequence in CSQ_CODING_HIGH_IMPACT:
                csq_order.append((loftee_filter, no_flag, consequence))

    # prioritization of mis and syn variant on protein coding transcripts
    csq_order.extend([(hl.null(hl.tstr), True, x) for x in
                      CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT])

    # Any variant on a non protein coding transcript (ie. where LOF = None)
    csq_order.extend([(hl.null(hl.tstr), True, x) for x in
                      CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT])

    csq_order = hl.literal({(x): i for i, x in enumerate(csq_order)})

    mt = mt.annotate_rows(**hl.sorted(mt.tx_annotation, key=lambda x: csq_order[
        (x.lof, hl.or_else(hl.is_missing(x.lof_flag), False), x.csq)])[0])

    return mt


def fix_loftee_beta_nonlofs(tc):
    keep_same = hl.literal(set(CSQ_CODING_HIGH_IMPACT))

    return tc.annotate(lof=hl.cond(keep_same.contains(tc.csq), tc.lof, hl.null('str')),
                       lof_flag=hl.cond(keep_same.contains(tc.csq), tc.lof_flag, hl.null('str')))


def get_baselevel_expression_for_genes(mt, gtex, gene_list, get_proportions = None,
                                       gene_maximums_kt_path = gtex_v7_gene_maximums_kt_path):

    gtex_table = gtex.key_by("transcript_id")
    genes = hl.literal(gene_list)

    # Filter context_ht to genes of interest
    mt = mt.annotate_rows(in_gene_of_interest=
                                     genes.find(lambda x: mt.vep.transcript_consequences.any(lambda tc: tc.gene_symbol == x)))
    mt = mt.filter_rows(mt.in_gene_of_interest != "NA")

    # Need to modify process consequences to ignore splice variants, because these can occur on intronic regions

    all_coding_minus_splice = list(set(all_coding_csqs) -
                                   set(['splice_acceptor_variant', 'splice_donor_variant','splice_region_variant' ]))


    def add_most_severe_consequence_to_consequence_minus_splice(
            tc: hl.expr.StructExpression) -> hl.expr.StructExpression:
        """
        Copied from gnomad_hail but slight change
        """

        csqs = hl.literal(all_coding_minus_splice)
        return tc.annotate(
            most_severe_consequence=csqs.find(lambda c: tc.consequence_terms.contains(c))
)
    # Add worst consequence within transcript consequences
    mt = (mt.annotate_rows(vep=mt.vep.annotate(
        transcript_consequences=mt.vep.transcript_consequences.map(add_most_severe_consequence_to_consequence_minus_splice))))

    # Explode on transcript consequences
    mt = mt.explode_rows(mt.vep.transcript_consequences)
    mt_kt = mt.rows()


    # Filter to positions in the CDS regions
    cds_intervals = hl.import_bed(
        "gs://gnomad-berylc/tx-annotation/hail2/browser_integration/gencode.v19.CDS.forHail2.bed")
    mt_kt = mt_kt.annotate(in_cds=hl.is_defined(cds_intervals[mt_kt.locus]))
    mt_kt = mt_kt.filter(mt_kt.in_cds)

    # Filter to protein coding transcripts only
    mt_kt = mt_kt.filter(mt_kt.vep.transcript_consequences.biotype == "protein_coding")

    # Filter to coding variants to only evalute those effects
    mt_kt = filter_table_to_csqs(mt_kt, all_coding_minus_splice)

    # To avoid double counting transcripts at a given base, key by transcript and position and dedup
    mt_kt = mt_kt.key_by(mt_kt.locus, mt_kt.vep.transcript_consequences.transcript_id)
    mt_kt = mt_kt.distinct()

    # Annotate mt with the gtex values (ie. join them)
    mt_kt = mt_kt.annotate(tx_data=gtex_table[mt_kt.vep.transcript_consequences.transcript_id])

    # Group by gene, symbol and position
    ht_sum_of_bases = mt_kt.group_by(ensg=mt_kt.vep.transcript_consequences.gene_id,
                                     symbol=mt_kt.vep.transcript_consequences.gene_symbol,
                                     locus=mt_kt.locus).aggregate(
        sum_per_base=hl.agg.array_sum(mt_kt.tx_data.agg_expression))

    tissue_ids = sorted([y.tissue.replace("-", "_").replace(" ", "_").replace("(", "_").replace(")", "_") for y in
                         gtex.values.take(1)[0]])
    d = {tiss: i for i, tiss in enumerate(tissue_ids)}

    ht_sum_of_bases = ht_sum_of_bases.annotate(
        **{tissue: ht_sum_of_bases.sum_per_base[d[tissue]] for tissue in tissue_ids})

    if get_proportions:
        gene_maximums_kt = hl.read_table(gene_maximums_kt_path)
        ht_sum_of_bases = ht_sum_of_bases.key_by(ht_sum_of_bases.locus)
        ht_sum_of_bases = ht_sum_of_bases.annotate(alleles = "fillter")
        ht_sum_of_bases = get_expression_proportion(tx_table = ht_sum_of_bases,
                                                    tissues_to_filter = ["sum_per_base"],
                                                    gene_maximum_kt  = gene_maximums_kt)
        ht_sum_of_bases = ht_sum_of_bases.key_by(ht_sum_of_bases.locus)
        ht_sum_of_bases = ht_sum_of_bases.drop(ht_sum_of_bases.alleles)

    return ht_sum_of_bases


