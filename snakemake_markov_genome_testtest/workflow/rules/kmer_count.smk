
rule kmer_count:
    input:
        f"results/{RUN_NAME}/{{kmer}}/markov_sim.fa",
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
    log:
        out=f"results/{RUN_NAME}/{{kmer}}/markov_kmc.log"
        err=f"results/{RUN_NAME}/{{kmer}}/markov_kmc.err"
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/markov_kmc.tsv"
    params: 
        k="{kmer}",
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp",
        dump=f"results/{RUN_NAME}/{{kmer}}/markov_kmc_dump",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"]
    shell:
        """
        mkdir -p {params.tmp}
        kmc -fm -k{params.k} -t{threads} {input} {params.dump} {params.tmp} 2>> {log.err} 1>> {log.out}
        kmc_tools transform {params.dump} dump {output.txt} 2>> {log.err} 1>> {log.out}
        rm -rf {params.tmp} {params.dump}*
        """


rule kmer_count_base:
    input:
        fasta_file=config["ref_genom"],
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
    log:
        out=f"results/{RUN_NAME}/{{kmer}}/base_kmc.log",
        err=f"results/{RUN_NAME}/{{kmer}}/base_kmc.err",
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/base_kmc.tsv",
    params: 
        k="{kmer}",
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp1",
        dump=f"results/{RUN_NAME}/{{kmer}}/base_dump",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"]
    shell:
        """
        mkdir -p {params.tmp}
        kmc -fm -k{params.k} -t{threads} {input.fasta_file} {params.dump} {params.tmp} 2>> {log.err} 1>> {log.out} ;
        kmc_tools transform {params.dump} dump {output.txt} 2>> {log.err} 1>> {log.out} ;
        rm -rf {params.tmp} {params.dump}*
        """


rule kmer_count_random:
    input:
        random=f"results/{RUN_NAME}/random/random_sim.fa",
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
    log:
        out=f"results/{RUN_NAME}/{{kmer}}/random_kmc.log",
        err=f"results/{RUN_NAME}/{{kmer}}/random_kmc.err",
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/random_kmc.tsv",
    params: 
        k="{kmer}",
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp2",
        dump=f"results/{RUN_NAME}/{{kmer}}/random_dump",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"]
    shell:
        """
        mkdir -p {params.tmp}
        kmc -fm -k{params.k} -t{threads} {input.random} {params.dump} {params.tmp} 2>> {log.err} 1>> {log.out}
        kmc_tools transform {params.dump} dump {output.txt} 2>> {log.err} 1>> {log.out}
        rm -rf {params.tmp} {params.dump}*
        """