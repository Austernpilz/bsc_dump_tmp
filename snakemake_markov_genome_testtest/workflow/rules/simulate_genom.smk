

rule generate_markov_genome:
    input:
        fasta_file=config["ref_genom"],
    params:
        k="{kmer}",
        s=config["bp"],
    output:
        f"results/{RUN_NAME}/{{kmer}}/markov_sim.fa",
    log:
        out=f"results/{RUN_NAME}/{{kmer}}/markov_sim.log",
        err=f"results/{RUN_NAME}/{{kmer}}/markov_sim.err",
    priority: 
        50
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/markov_sim.tsv",
    message:
        """--- Generating genome sequence from {wildcards.kmer}-kmer""",
    shell:
        """
        markov_genome \
        --input {input.fasta_file} \
        --output {output} \
        --order {params.k} \
        --lens {params.s} 2>> {log.err} 1>> {log.out} ;
        /scripts/fix_broken_files.sh {output} 
        """

