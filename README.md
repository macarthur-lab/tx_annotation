# Transcript expression-aware annotation 

Welcome to the repository for the "Transcript expression-aware annotation improves rare variant discovery and interpretation" manuscript. Here we'll outline how to get transcript expression values for your variant file and isoform expression expression matrix of interest, and outline the commands and code to recreate analyses in the pre-print. 

## Applying transcript expression aware annotation to your own dataset

You will need 
  1) A variant file that has columns for chrom, pos, ref and alt
  2) An isoform expression matrix 
  3) The ability to use Hail locally or on a cloud platform 
  
You can have additional columns in your variant file, which will be maintained, and only new columns of transcript-expression annotation will be added. Your isoform expression matrix must start with two columns : 1. transcript_id and 2. gene_id. The remaining columns can be any samples or tissues. If you have biological replicates, they should be numbered with a '.' delimiter (e.g. MuscleSkeletal.1, MuscleSkeletal.2, MuscleSkeletal.3). 
 
Instructions to set up Hail [can be found in the Hail docs](https://hail.is/docs/0.2/)

If you're unable to set up Hail in your local environment, we have released the pext values for every possible SNV in the genome in both a tsv and HT format in the folder: gs://gnomad-public/papers/2019-tx-annotation/pre_computed/

More information about this file's format is below.

Please be aware that while we don't expect any issues, the files may be iterated upon until publication of the manuscript! 

We will walk through an example of annotating *de novo* variants in autism and developmental delay / intellectual disability with the GTEx v7 dataset. 

#### 0) Start a cluster and a Hail environment
We recommend using [hailctl dataproc from the Hail team](https://hail.is/docs/0.2/hail_on_the_cloud.html) for Google Cloud. Note that this is a successor to cloudtools from Ben Neale's lab.

You will need the gnomAD and tx-annotation init scripts, which are both publically available. To start a cluster:

```
hailctl dataproc start tutorial --init gs://gnomad-public/tools/inits/master-init.sh,gs://gnomad-public/papers/2019-tx-annotation/tx-annotation-init.sh --num-preemptible-workers 8 --vep GRCh37
```

To connect to the clustter 
```
hailctl dataproc connect tutorial nb
```
And start a Hail Jupyter notebook

#### 1) Prepare the variant file 
The variant file we'll be using for the tutorial is available at : gs://gnomad-public/papers/2019-tx-annotation/data/de_novo_variants/asd_ddid_de_novos.txt

This is what the first line of the file looks like : 
> DataSet CHROM POSITION  REF ALT	GENE_NAME	VEP_functional_class_canonical	MPC	loftee	group
> ASC_v15_VCF 1 94049574  C A BCAR3 splice_donor_variant  NA	HC	Control

Again, we will only use the chrom, pos, ref, alt columns, and will add additional columns. The VEP columns in the file are based on the canonical trasncript, so we will re-VEP. 

Note that unless the variant file you're using is already available as HT that's been vep'd in Hail, we recommend running VEP since the pext code uses the nested VEP format for calculation. 

This is how to import the file into Hail, define the variant field, vep, and write the MT. Note that this VEP configuration will also annotate with LOFTEE v.1.0

0 - At the top of your script specify `from tx_annotation import *` which will start a Hail environment, and import necessary parts of the gnomAD repository.

1 - Import file as a table 
```python
rt = hl.import_table("gs://gnomad-public/papers/2019-tx-annotation/data/de_novo_variants/asd_ddid_de_novos.txt")
```
2 - Define the variant in terms of chrom:pos:ref:alt and have Hail parse it, which will create locus and alleles fields
```python
rt = rt.annotate(variant=rt.CHROM + ':' + rt.POSITION + ":" + rt.REF + ":" + rt.ALT)
rt = rt.annotate(** hl.parse_variant(rt.variant))
rt = rt.key_by(rt.locus, rt.alleles)
```
3 - Make a MT from the Table, and repartition for speed (rule of thumb is ~2k variants per partition)
```python
mt = hl.MatrixTable.from_rows_table(rt)
mt = mt.repartition(10)
```
4 - VEP and write out the MT
```python
annotated_mt = hl.vep(mt, vep_config)
annotated_mt.write("gs://gnomad-public/papers/2019-tx-annotation/results/de_novo_variant/asd_ddid_de_novos.vepped.021819.mt")
```

#### 2) Prepare the isoform expression file 

We'll use the GTEx v7 isoform expression file. Here is what the header of the file looks like 
> transcript_id   gene_id GTEX-1117F-0226-SM-5GZZ7        GTEX-1117F-0426-SM-5EGHI        GTEX-1117F-0526-SM-5EGHJ        
> ENST00000373020.4       ENSG00000000003.10      26.84   4.13    13.54  

We've replaced the sample names with unique tissue names, so that samples with the same tissue are labelled as WholeBlood.1, WholeBlood.2, WholeBlood.3 etc:
> transcript_id   gene_id Adipose-Subcutaneous.1  Muscle-Skeletal.2       Artery-Tibial.3 
> ENST00000373020.8       ENSG00000000003.14      26.32   3.95    13.23   

We first need to get the median expression of all transcripts per tissue. This can be carried out using the `get_gtex_summary()` function (the function name is a misnomer, as it can work on non-GTEx files).

```python
gtex_isoform_expression_file = /path/to/text/file/with/isoform/quantifications
gtex_median_isoform_expression_mt = /path/to/matrix_table/file/you/want/to/create
get_gtex_summary(gtex_isoform_expression_file,gtex_median_isoform_expression_mt )

```
If you'd like to get mean isoform expression accross tissues and not median, add get_medians = False to the command. If you want to also export the median isoform expression per tissue file as a tsv, add make_per_tissue_file = True 

Unfortunately, we can't share the per-sample GTEx RSEM file as it requires dbGAP approval. However, running this on the GTEx v7 dataset creates: gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.V7.tx_medians.110818.mt which is the file used for the analyses in the manuscript and the file you can use for your annotation GTEx v7 annotation. 

At this point, you'll also need a separate file with gene expression values per tissue, with the tissue names matching the median isoform expression file. For the manuscript, we directly imported gene expression values provided by GTEx, which were created using RNASeQC, from the GTEx portal website. They are available here: gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.v7.gene_expression_per_gene_per_tissue.120518.kt

#### 3) Add pext values
All you have to do at this point is import your VEP'd variant matrix table, and run the tx_annotate() function!

1 - Import VEP'd variant MT, and median isoform expression MT: 
```python
mt, gtex = read_tx_annotation_tables(ddid_asd_de_novos, gtex_v7_tx_summary_mt_path, "mt")
```
2 - Run tx_annotation
```python
ddid_asd = tx_annotate_mt(mt, gtex,
                          tx_annotation_type = "proportion",
                          filter_to_csqs=all_coding_csqs)

```

This command by default will remove certain GTEx tissues with <100 samples, reproductive tissues, or cell lines (specified in `tx_annotation_resources` and in the manuscript). 

- If you don't want to remove these tissues (or if you are not working with GTEx) specify `tissues_to_filter = None`.
- If you'd like to get the non-normalized ext values instead of pext, specify `tx_annotation_type = "expression"`.
- Not specifying `filter_to_csqs=all_coding_csqs` will add pext values to  non-coding variants (which may be desired behavior based on your goals). Note that splice variants are considered coding variants here. The description of coding csqs is available in `tx_annotation_resources` as `all_coding_csqs`.
- If you're only interested in getting pext for a certain group of genes, you can specify that with `filter_to_genes`. This will return the same file, but will only add the pext values to the genes of interest. An example of adding pext while specifying genes is below (under the ClinVar - gnomAD comparison section) 

The function returns your variant MT with a new field called `tx_annotation` (ddid_asd_de_novos_with_pext above). 

At this point, you can choose what annotation you want to use for a given variant (for example, you may be interested in any pLoF variant, or variants found on certain set of transcripts, or just variants found on the canonical transcript - the last of which  sort of defeats the point of using this method). In the manuscript we used the worst consequence accross transcripts, which is the context for which we see this method being most powerful. If you'd also like to use the worst consequence, and pull out pext values for the worst consequence, we have helper functions available: 

#### 4) Optional post-processing to pull out pext values for the worst consequence annotation

At this point you will remove all variants that did not receive a pext value (e.g. if you specific `filter_to_csqs = all_coding_csqs` this will remove noncoding variants). At this point, we don't support the OS annotation in LOFTEE, which add pLoF annotations to missense and synonymous variants (for example, a synonymous variant can be called LOFTEE HC in the latest LOFTEE release if it's predicted to affect splicing). We therefore replace OS annotations with the original annotation (ie. we replace the HC for a synonymous variant with ""). Finally, we extract the worst consequence, and create one column per tissue. 

1 - Remove variants that did not receive a pext annotation (ie. noncoding variants)
```python
ddid_asd = ddid_asd.filter_rows(~hl.is_missing(ddid_asd.tx_annotation))

```
2 - Overwrite LOFTEE OS variants with original variant annotation 
```python
ddid_asd = ddid_asd.annotate_rows(tx_annotation=ddid_asd.tx_annotation.map(fix_loftee_beta_nonlofs))
```
3 - Pull out worst consequence
```python
ddid_asd = pull_out_worst_from_tx_annotate(ddid_asd)
```

At this point you can write out the file with `ddid_asd.rows().export("out_file")`

This will create the transcript annotated *de novo* variant file used in Figure 4 of the manuscript. We've exported the result of this code snippet here: 
gs://gnomad-public/papers/2019-tx-annotation/results/de_novo_variant/asd_ddid_de_novos.tx_annotated.proprotion.021819.tsv.bgz


## Analyses in manuscript 
Here we'll detail the commands for obtaining pext values for some of the analyses in manuscript. This will go over the analysis of:
- Getting baselevel expression values for a gene (Figure 2B) 
- Comparison of highly conserved and unconserved regions (Figure 3A) 
- Comparison of % variant filtered with pext < 0.1 in haploinsufficient disease genes (Figure 4A)

Note that scripts for these and other analyses are availabel in `/analyses/` folder in this repository. The paths to the files are available in `tx_annotation_resources.py` If you find something is missing, please e-mail me at berylc@broadinstitute.org


#### Getting baselevel expression values 
The idea here is that you annotate the expression of a given *position* as opposed to a variant consequence pair. The baselevel pext value will always be higher than any of the variant annotation pext values, because the base value is just the sum of the expression of protein coding transcripts that overlap the coding base. This baselevel value is what we show in the gnomAD browser. Just because a position has a high baselevel value though, *does not* mean that say, a pLoF at that position has a high pLoF value.

We get these baselevel values by using the sites table of all possible variant in the genome. We sum of the expression of all transcripts overlapping that base, where there's a coding consequence. 

- TCF4 
```python
from tx_annotation import * 
mt, gtex = read_tx_annotation_tables(context_ht_path, gtex_v7_tx_summary_mt_path, "ht")
gene_baselevel= get_baselevel_expression_for_genes(mt, gtex, gene_list = {'TCF4'})
gene_baselevel.export("gs://gnomad-public/papers/2019-tx-annotation/results/TCF4.baselevel.ext.021319.tsv.bgz")
```

The resulting file is used in Fig2B and Supp Fig 4. 

You can specify any number of genes you want in `gene_list`. If you don't specify any genes, it will annotate all positions in the exome. 

- SCN2A using fetal isoform expression 
```python
hbdr_fetal_path = "gs://gnomad-public/papers/2019-tx-annotation/data/HBDR.RSEM.sample_specific.tx_medians.021719.mt"
mt, hbdr_fetal = read_tx_annotation_tables(context_ht_path, hbdr_fetal_path, "ht")
gene_baselevel= get_baselevel_expression_for_genes(mt, hbdr_fetal, gene_list = {'SCN2A'})

```
This file was used in Supp Fig 6D.

#### Comparison of pext in highly conserved and unconserved regions


1 - Read in baselevel expresison and phyloCSF files
```python
phylocsf = hl.import_table(phylocsf_file_path, impute = True)
all_baselevel_ht = hl.read_table(all_baselevel_ht_path)
phylocsf = phylocsf.annotate(chrom = phylocsf.chromosome_name.replace("chr",""))
```
Note that all_baselevel_ht_path is the file used to create the tx_annotation tracks in the [gnomAD browser](https://gnomad.broadinstitute.org/)

2 - Define regions of high and low conservation, and filter remaning regions
```python
phylocsf = phylocsf.annotate(conservation_type = hl.case(missing_false=True)
                             .when(phylocsf.max_score > 1000, "high")
                             .when(phylocsf.max_score < -100, "low")
                             .default('filter'))
phylocsf = phylocsf.filter(phylocsf.conservation_type != "filter")
```

3 -  Make intervals in the phyloCSF file
phylocsf = phylocsf.annotate(chrom = phylocsf.chromosome_name.replace("chr",""))
```python
phylocsf = phylocsf.annotate(
    interval = hl.interval(hl.locus(phylocsf.chrom, phylocsf.start_coordinate), hl.locus(phylocsf.chrom, phylocsf.end_coordinate)),
    interval_name = phylocsf.chrom + ":" + hl.str(phylocsf.start_coordinate) + "-" + hl.str(phylocsf.end_coordinate) )
phylocsf = phylocsf.key_by(phylocsf.interval)
```
4 - Filter the baselevel expression file to the intervals of high or low conservation in the phyloCSF file
```python
all_baselevel_ht = all_baselevel_ht.annotate(**phylocsf[all_baselevel_ht.locus])
all_baselevel_ht= all_baselevel_ht.filter(hl.is_defined(all_baselevel_ht.conservation_type), keep=True)
```
 5- Get mean pext in these intervals
``` python
mean_proportion_in_interval = (all_baselevel_ht.group_by(symbol = all_baselevel_ht.symbol,
                                                         ensg = all_baselevel_ht.ensg,
                                                         enst = all_baselevel_ht.transcript_id,
                                                         interval = all_baselevel_ht.interval_name,
                                                         conservation_type = all_baselevel_ht.conservation_type).
                               aggregate(mean_of_mean_pext =
                                        hl.agg.filter(~hl.is_nan(all_baselevel_ht.mean_prop_conservation),
                                                      hl.agg.mean(all_baselevel_ht.mean_prop_conservation))))
```
6 - Export the file for plotting
```python
mean_proportion_in_interval.export("gs://gnomad-public/papers/2019-tx-annotation/results/conservation.phylocsf.vs.pext.021219.tsv.bgz")
```

#### Comparison of % variant filtered with pext < 0.1 in haploinsufficient disease genes (Figure 4A)

This also serves as an example of annotating only a subset of genes in a variant table. 

Here we will annotate variants in HI genes in the gnomAD exomes sites HT, the gnomAD genomes sites HT, and the ClinVar HT with pext values. 

```
out_dir = "gs://gnomad-public/papers/2019-tx-annotation/results/gene_list_comparisons/"
```

1 - Import HI genes 
```python
hi_genes = import_gene_list(curated_haploinsufficient_genes, gene_column="ENSGID", ensg=True)
```
There are two options for importing gene lists, either importing ENSG IDs, or importing gene symbols. `gene_column` refers to the column in the file that contains your gene names, if the values are ENSGs, specify `ensg = True`. You can specify `peek = True` if you'd like to just like to import the gene list file and take a peek without doing anything.

2 - Annotate gnomAD exomes
```python
mt, gtex = read_tx_annotation_tables(gnomad_release_mt_path, gtex_v7_tx_summary_mt_path, "ht")
mt = mt.filter_rows(hl.len(mt.filters) == 0)
mt_gnomad_hi = tx_annotate_mt(mt, gtex,"proportion",
                              filter_to_csqs=lof_csqs,
                              filter_to_genes=hi_genes, gene_column_in_mt="gene_id")
mt_gnomad_hi = mt_gnomad_hi.filter_rows(~hl.is_missing(mt_gnomad_hi.tx_annotation))
mt_gnomad_hi = pull_out_worst_from_tx_annotate(mt_gnomad_hi)
mt_gnomad_hi.rows().export("%sHI_genes.gnomad.exomes.r2.1.tx_annotated.021519.tsv.bgz" %out_dir)

- `gene_column_in_mt` is one of either `gene_id` (ENSG) or `gene_symbol` and tells the function which VEP field to look to filter to genes.

- `mt = mt.filter_rows(hl.len(mt.filters) == 0)` filters variants to only those that are RF PASS.

```
3 - Annotate gnomAD genomes 
```python
mt_genomes, gtex = read_tx_annotation_tables(gnomad_genomes_release_mt_path, gtex_v7_tx_summary_mt_path, "ht")
mt_genomes = mt_genomes.filter_rows(hl.len(mt_genomes.filters) == 0)
mt_gnomad_genomes_hi = tx_annotate_mt(mt_genomes, gtex,"proportion",
                              filter_to_csqs=lof_csqs,
                              filter_to_genes=hi_genes, gene_column_in_mt="gene_id")
mt_gnomad_genomes_hi = mt_gnomad_genomes_hi.filter_rows(~hl.is_missing(mt_gnomad_genomes_hi.tx_annotation))
mt_gnomad_genomes_hi = pull_out_worst_from_tx_annotate(mt_gnomad_genomes_hi)
mt_gnomad_genomes_hi.rows().export("%sHI_genes.gnomad.genomes.r2.1.tx_annotated.021619.tsv.bgz" %out_dir)
```

4 - Annotate ClinVar 
```python
clinvar_mt, gtex = read_tx_annotation_tables(clinvar_ht_path, gtex_v7_tx_summary_mt_path, "ht")
mt_clinvar_hi = tx_annotate_mt(clinvar_mt, gtex,"proportion",
                               filter_to_csqs=lof_csqs, filter_to_genes=hi_genes,
                               gene_column_in_mt="gene_id")
mt_clinvar_hi = mt_clinvar_hi.filter_rows(~hl.is_missing(mt_clinvar_hi.tx_annotation))
mt_clinvar_hi = pull_out_worst_from_tx_annotate(mt_clinvar_hi)
mt_clinvar_hi = mt_clinvar_hi.annotate_rows(**mt_clinvar_hi.info)
mt_clinvar_hi = mt_clinvar_hi.drop("vep", "tx_annotation","info") 
mt_clinvar_hi.rows().export("%sHI_genes.clinvar.alleles.single.b37.tx_annotated.021519.tsv.bgz" %out_dir)
```
The fields are dropped to save space. 
