# Looking for Sdic transcripts in ISO-Seq data

Jose Ranz asked us to look in our testis Iso-Seq data for potential Sdic isoforms.

> ... we would like to have a better, more precise information about the gene
> models of the different copies.  We have thought that Iso-seq data from
> testes or whole-bodies would be ideal ...

> we would like to have is all transcripts that carry a particular exon...
> Would it be possible to ... blast search and pull out all the transcripts in
> a multifasta format?  

Jose provided this query string (10/10/2019):

GATTGGACAATGGGCTTAGTACTGATTAAGTTTTTACGATCAACGTATTCTACTTTG

# Outputs

We have Iso-Seq data for:

* Adult whole body (including gonads)
* Adult testis
* Third instar larval testis

Outputs include:

* PSL file showing alignment of query to various full length transcripts
* FASTQ file of full length transcripts that matched the Sdic query sequence

# To repeat analysis

First you need to obtain the data (polished FASTA and FASTQ) and fix the directory structure in Snakefile.

## Create MiniConda environment

```bash
conda env create --file conda_environment.yaml -p ./env
conda activate ./env
```

## Re-run workflow

```bash
snakemake -p
```
