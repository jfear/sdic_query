from snakemake.io import expand

SDIC_EXON = "GATTGGACAATGGGCTTAGTACTGATTAAGTTTTTACGATCAACGTATTCTACTTTG"

patterns = {
    "adult_testis": "/home/fearjm/local_data_store/dmel_pacbio/output/pacbio-wf/w1118_testi1/w1118_testi1.polished.hq.fasta.gz",
    "adult_whole_body": "/home/fearjm/local_data_store/dmel_pacbio/output/pacbio-wf/w1118_wmal1/w1118_wmal1.polished.hq.fasta.gz",
    "larval_testis": "/home/fearjm/local_data_store/dmel_larval_pacbio/output/l3_testes/l3_testes.polished.hq.fasta.gz",
}


rule all:
    input: expand("output/{sample}_sdic_reads.fastq.gz", sample=["adult_testis", "adult_whole_body", "larval_testis"])


rule query_to_fasta:
    """Convert SDIC Exon Seq to FASTA"""
    output: temp("output/query.fa")
    params: query = SDIC_EXON
    shell: """
        echo '>sdic query' > {output[0]} \
        && echo {params.query} >> {output[0]}
    """


def get_fastq_name(wildcards):
    return patterns[wildcards.sample]


rule blat:
    """Use blat to find Full length transcripts with sequence."""
    input: 
        ref = get_fastq_name,
        query = rules.query_to_fasta.output[0]
    output: "output/{sample}.psl"
    shell: "blat {input.ref} {input.query} {output[0]}"


rule unzip_fastq:
    input: lambda wildcards: get_fastq_name(wildcards).replace("fasta", "fastq")
    output: temp("output/{sample}.fastq")
    shell: "gunzip -c {input} > {output[0]}"


rule pull_out_sdic_reads:
    """Use blat results to pull out full length transcript sequences"""
    input:
        fq = "output/{sample}.fastq",
        psl = rules.blat.output[0]
    output: temp("output/{sample}_sdic_reads.fastq")
    run:
        import pandas as pd
        from Bio import SeqIO
        transcript_ids = (
            pd.read_csv(input.psl, skiprows=5, sep="\t", header=None)
            .iloc[:, 13]
            .values.tolist()
        )
        sdic_records = [
            record
            for record in SeqIO.parse(input.fq, "fastq")
            if record.id in transcript_ids
        ]

        with open(output[0], "w") as fh:
            SeqIO.write(sdic_records, fh, "fastq")


rule zip_sdic_reads:
    input: rules.pull_out_sdic_reads.output[0]
    output: "output/{sample}_sdic_reads.fastq.gz"
    shell: "gzip {input}"
