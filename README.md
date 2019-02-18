# Transcript expression-aware annotation 

Welcome to the repository for the "Transcript expression-aware annotation improves rare variant discovery and interpretation" manuscript. Here we'll outline how to get transcript expression values for your variant file and isoform expression expression matrix of interest, and outline the commands and code to recreate analyses in the pre-print. 

## Applying transcript expression aware annotation to your own dataset

You will need 
  1) A variant file that has columns for chrom, pos, ref and alt
  2) An isoform expression matrix 
  3) The ability to use Hail locally or on a cloud platform 
  
You can have additional columns in your variant file, which will be maintained, and only new columns of transcript-expression annotation will be added. Your isoform expression matrix must start with two columns : 1. transcript id and 2. gene_id. The remaining columns can be any samples or tissues. If you have biological replicates, they should be numbered with a '.' delimiter (e.g. MuscleSkeletal.1, MuscleSkeletal.2, MuscleSkeletal.3). 
 
Instructions to set up Hail are laid out in the [Hail docs](https://hail.is/docs/0.2/)

If you're unable to set up Hail in your local environment, we have also release the pext values for every possible SNV in the genome: gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021819.tsv.bgz

More information about this files format is below.

Please be aware that while we don't expect any issues, the files may be iterated upon until publication of the manuscript! 

We will walk through an example of annotating *de novo* variants in autism and developmental delay / intellectual disability with the GTEx v7 dataset. 

#### 0) Start a cluster and a Hail environment
We recommend using [cloud tools from Neale lab](https://github.com/Nealelab/cloudtools) for Google Cloud.

You will need the gnomAD and tx-annotation init scripts, which are both publically available. Then you can start a cluster:

```
cluster start tutorial --worker-machine-type n1-highmem-8 --spark 2.2.0 --version 0.2 --init gs://gnomad-public/tools/inits/master-init.sh,gs://gnomad-public/papers/2019-tx-annotation/tx-annotation-init.sh --num-preemptible-workers 8
```

At the top of your script specify `from tx_annotation import *` which will start a Hail environment, and import necessary parts of the gnomAD repository.

#### 1) Prepare the variant file 
The variant file is available at : gs://gnomad-public/papers/2019-tx-annotation/data/asd_ddid_de_novos.txt

This is what the first line of the file looks like : 
>      DataSet CHROM POSITION REF ALT Child_ID Child_Sex GENE_NAME VEP_functional_class_canonical MPC loftee   group
>       ASC_v15_VCF     1 94049574   C   A 13069.s1      Male     BCAR3           splice_donor_variant  NA     HC Control

Again, we will only use the chrom, pos, ref, alt columns, and will add additional columns. 

In order to add pext values, you must annotate with VEP. This is how to import the file into Hail, define the variant field, vep, and write the MT. 

1 - Import file as a table 
```
rt = hl.import_table("gs://gnomad-public/papers/2019-tx-annotation/data/asd_ddid_de_novos.txt")
```
2 - Define the variant in terms of chrom:pos:ref:alt and have Hail parse it, which will create a locus and alleles column
```
rt = rt.annotate(variant=rt.CHROM + ':' + rt.POSITION + ":" + rt.REF + ":" + rt.ALT)
rt = rt.annotate(** hl.parse_variant(rt.variant))
rt = rt.key_by(rt.locus, rt.alleles)
```
3 - Make a MT from table, and repartition for speed (rule of thumb is ~2k variants per partition)
```
mt = hl.MatrixTable.from_rows_table(rt)
mt = mt.repartition(10)
```
4 - VEP and write out the MT
```
annotated_mt = hl.vep(mt, vep_config)
annotated_mt.write("gs://gnomad-public/papers/2019-tx-annotation/results/de_novo_variants/asd_ddid_de_novos.vepped.021819.mt")
```

#### 2) Prepare the isoform expression file 

We'll use the GTEx v7 isoform expression file. Here is what the header of the file looks like 
> transcript_id   gene_id GTEX-1117F-0226-SM-5GZZ7        GTEX-1117F-0426-SM-5EGHI        GTEX-1117F-0526-SM-5EGHJ        
> ENST00000373020.4       ENSG00000000003.10      26.84   4.13    13.54  

We've replaced the sample names with unique tissue names, so that samples with the same tissue are labelled as WholeBlood.1, WholeBlood.2, WholeBlood.3 etc:
> transcript_id   gene_id Adipose-Subcutaneous.1  Muscle-Skeletal.2       Artery-Tibial.3 
> ENST00000373020.8       ENSG00000000003.14      26.32   3.95    13.23   

We first to get the median expression of all transcripts per tissue. This can be carried out using the `get_gtex_summary()` function (the function name is a bit of a misnomer, as it can work on non-GTEx files).

```
gtex_isoform_expression_file = /path/to/text/file/with/isoform/quantifications
gtex_median_isoform_expression_mt = /path/to/matrix_table/file/you/want/to/create
get_gtex_summary(gtex_isoform_expression_file,gtex_median_isoform_expression_mt )

```
If you'd like to get mean isoform expression accross tissues and not mean, add get_medians = False to the command. If you want to also export the median isoform expression per tissue file as a tsv, add make_per_tissue_file = True 

Unfortunately, we can't share the per-sample GTEx RSEM file as it requires dbGAP approval. However, running this on the GTEx v7 dataset creates: gs://gnomad-public/papers/2019-tx-annotation/data/GTEx.V7.tx_medians.110818.mt which is the file used for the analyses in the manuscript. 

#### 3) Add pext values


## Analyses in manuscript 

###### Conservation and pext
d
