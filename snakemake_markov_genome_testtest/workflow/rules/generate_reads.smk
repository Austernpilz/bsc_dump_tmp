
rule simlord_random:
    output:
        fasta=f"results/{RUN_NAME}/random/random_sim.fa",
        reads=temp(f"results/{RUN_NAME}/random/sample.fastq")
    params:
        s=config["bp"],
        rand=f"results/{RUN_NAME}/random/",
        reads=temp(f"results/{RUN_NAME}/random/sample")
    benchmark:
        f"results/{RUN_NAME}/random/random_sim.tsv",
    message:
        """--- generating random reads with simlord."""
    shell:
        """
        mkdir -p {params.rand}
        simlord --generate-reference 0.5 {params.s} --save-reference {output.fasta} \
        -c 0.9 --no-sam {params.reads} ;
        /scripts/fix_broken_files.sh {output.reads} ;
        /scripts/fix_broken_files.sh {output.fasta} ;
        """


rule simlord_reads_baseline:
    input:
        fasta_file=config["ref_genom"],
    output:
        reads=temp(f"results/{RUN_NAME}/base/sample.fastq")
    params:
        reads=temp(f"results/{RUN_NAME}/base/sample")
    message:
        """--- generating baseline reads with simlord."""
    shell:
        """
        simlord --read-reference {input.fasta_file} -c 0.9 --no-sam {params.reads}
        /scripts/fix_broken_files.sh {output.reads}
        """

rule simlord_reads:
    input:
        f"results/{RUN_NAME}/{{kmer}}/markov_sim.fa"
    output:
        reads=temp(f"results/{RUN_NAME}/{{kmer}}/sample.fastq")
    params:
        reads=temp(f"results/{RUN_NAME}/{{kmer}}/sample")
    message:
        """--- generating reads with simlord."""
    shell:
        """
        simlord --read-reference {input} -c 0.9 --no-sam {params.reads}
        /scripts/fix_broken_files.sh {output.reads}
        """

