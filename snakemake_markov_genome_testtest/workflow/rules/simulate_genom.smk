

rule generate_markov_genome:
    input:
        fasta_file=config["ref_genom"],
    params:
        k={kmer},
        s=config["bp"],
    output:
        temp(f"results/{RUN_NAME}/{{kmer}}/markov_sim.fa"),
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/markov_sim.tsv",
    message:
        """--- Generating genome sequence from {wildcards.kmer}-kmer""",
    resources: config["mb"] * 4
    priority: 50
    shell:
        """
        markov_genome \
        --input {input.fasta_file} \
        --output {output} \
        --order {params.k} \
        --lens {params.s} 2>> {log.err} 1>> {log.out} ;
        """

