

rule kmer_count_random_1234:
    input:
        wone=f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
        random=f"results/{RUN_NAME}/random/random_sim.fa",
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/random_kmc.tsv",
    params:
        k="{kmer}",
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp/random/",
        dump=f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/random_dump",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"] * confi["kmc_t"] + 10000
    shell:
        """
        kmc -fm -k{params.k} -t{threads} {input.random} {params.dump} {params.tmp}
        kmc_tools transform {params.dump} dump {output.txt}
        """