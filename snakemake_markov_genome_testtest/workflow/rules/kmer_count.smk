
rule look_out_for_folders:
    input:
        f"results/{RUN_NAME}/{{kmer}}/markov_sim.fa",
        fasta_file=config["ref_genom"],
        random=f"results/{RUN_NAME}/random/random_sim.fa",
    output:
        f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
    priority:
        1000
    params:
        all_mer = [ 
            f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/", f"results/{RUN_NAME}/{{kmer}}/tmp/markov/", 
            f"results/{RUN_NAME}/{{kmer}}/tmp/random/", f"results/{RUN_NAME}/{{kmer}}/tmp/base/",
            f"results/{RUN_NAME}/{{kmer}}/tmp/union_markov/", f"results/{RUN_NAME}/{{kmer}}/tmp/union_random/",
            f"results/{RUN_NAME}/{{kmer}}/tmp/inter_markov/", f"results/{RUN_NAME}/{{kmer}}/tmp/inter_random/",
        ]
    run:
        for p in all_mer:
            if not os.path.isdir(p):
                os.makedirs(p, exist_ok=True)
        with open(output, 'w+') as f:
            f.write("done")
        print("dirs are created")


rule kmer_count_markov:
    input:
        f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
        fasta = f"results/{RUN_NAME}/{{kmer}}/markov_sim.fa",
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
        dump=[f"results/{RUN_NAME}/{{kmer}}/markov_kmc_dump.kmc_suf", f"results/{RUN_NAME}/{{kmer}}/markov_kmc_dump.kmc_pre"]
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/markov_kmc.tsv"
    params: 
        k="{kmer}",
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp/markov/",
        dump=f"results/{RUN_NAME}/{{kmer}}/markov_kmc_dump",
    threads: config["kmc_t"],
    resources: mem_mb=config["mb"] * config["kmc_t"] + 10000
    shell:
        """
        kmc -fm -k{params.k} -t{threads} {input.fasta} {params.dump} {params.tmp}
        kmc_tools transform {params.dump} dump {output}
        """


rule kmer_count_base:
    input:
        f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
        fasta_file=config["ref_genom"],
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
        dump=[f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump.kmc_pre",f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump.kmc_suf"]
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/base_kmc.tsv",
    params: 
        k="{kmer}",
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp/base/",
        dump=f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"] * config["kmc_t"] + 10000
    shell:
        """
        kmc -fm -cs104875 -k{params.k} -t{threads} {input.fasta_file} {params.dump}_dump {params.tmp}
        kmc_tools transform {params.dump}_dump dump {output.txt}
        """


rule kmer_count_random:
    input:
        f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
        random=f"results/{RUN_NAME}/random/random_sim.fa",
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
        dump=[f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/random_dump.kmc_pre",f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/random_dump.kmc_surf"]
    benchmark:
        f"results/{RUN_NAME}/{{kmer}}/random_kmc.tsv",
    params: 
        k="{kmer}",
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp/random/",
        dump=f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/random_dump",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"] * config["kmc_t"] + 10000
    shell:
        """
        kmc -fm -k{params.k} -t{threads} {input.random} {params.dump} {params.tmp}
        kmc_tools transform {params.dump} dump {output.txt}
        """



rule kmer_count_intersect_markov:
    input:
        fasta = f"results/{RUN_NAME}/{{kmer}}/markov_sim.fa",
        fasta_file=config["ref_genom"],
        f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
        dump_markov=[f"results/{RUN_NAME}/{{kmer}}/markov_kmc_dump.kmc_suf", f"results/{RUN_NAME}/{{kmer}}/markov_kmc_dump.kmc_pre"]
        dump_base=[f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump.kmc_pre",f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump.kmc_suf"]
    output:
        dump_inter = [f"results/{RUN_NAME}/{{kmer}}/tmp/inter_markov.kmc_pre",f"results/{RUN_NAME}/{{kmer}}/tmp/inter_markov.kmc_suf",]
        txt=f"results/{RUN_NAME}/{{kmer}}/intersect_markov_base.txt",
    params: 
        k="{kmer}",
        tmp = f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/",
        dump0 = f"results/{RUN_NAME}/{{kmer}}/markov_kmc_dump",
        dump1 = f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump",
        dump2 = f"results/{RUN_NAME}/{{kmer}}/tmp/inter_markov",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"] * config["kmc_t"] + 10000
    shell:
        """
        kmc_tools simple {params.dump0} {params.dump1} {params.dump2}
        kmc_tools transform {params.dump2} dump {output.txt}
        """
# kmc -fm -k{params.k} -t{threads} {input.fasta} {params.dump0} {params.tmp}
#         kmc -fm -k{params.k} -t{threads} {input.fasta_file} {params.dump1} {params.tmp}

rule kmer_count_intersect_random:
    input:
        random=f"results/{RUN_NAME}/random/random_sim.fa",
        fasta_file=config["ref_genom"],
        f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
        dump_base=[f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump.kmc_pre",f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump.kmc_suf"]
        dump_random=[f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/random_dump.kmc_pre",f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/random_dump.kmc_surf"]
    output:
        txt=f"results/{RUN_NAME}/{{kmer}}/intersect_random_base.txt",
        dump_inter = [f"results/{RUN_NAME}/{{kmer}}/tmp/union_markov/inter.kmc_suf",f"results/{RUN_NAME}/{{kmer}}/tmp/union_markov/inter.kmc_pre",]
    params: 
        k="{kmer}",
        tmp = f"results/{RUN_NAME}/{{kmer}}/tmp/random/",
        dump0 = f"results/{RUN_NAME}/{{kmer}}/tmp/tmp/random_dump",
        dump1 = f"results/{RUN_NAME}/{{kmer}}/tmp/base_dump",
        dump2 = f"results/{RUN_NAME}/{{kmer}}/tmp/union_markov/inter",
    threads: config["kmc_t"]
    resources: mem_mb=config["mb"] * config["kmc_t"] + 10000
    shell:
        """
        kmc -fm -k{params.k} -t{threads} {input.random} {params.dump0} {params.tmp}
        kmc -fm -k{params.k} -t{threads} {input.fasta_file} {params.dump1} {params.tmp}
        kmc_tools simple {params.dump0} {params.dump1} {params.dump2}
        kmc_tools transform {params.dump2} dump {output.txt}
        """

rule close_tmp_again:
    input:
        dump_inter0 = [f"results/{RUN_NAME}/{{kmer}}/tmp/inter_markov.kmc_pre",f"results/{RUN_NAME}/{{kmer}}/tmp/inter_markov.kmc_suf",]
        dump_inter1 = [f"results/{RUN_NAME}/{{kmer}}/tmp/union_markov/inter.kmc_suf",f"results/{RUN_NAME}/{{kmer}}/tmp/union_markov/inter.kmc_pre",]
        f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
        f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/intersect_markov_base.txt",
        f"results/{RUN_NAME}/{{kmer}}/intersect_random_base.txt",
        f"results/plots/kmer_{RUN_NAME}_{{kmer}}.tsv",
        f"results/plots/kmer_{RUN_NAME}_{{kmer}}_intersect.tsv",
    params:
        tmp=f"results/{RUN_NAME}/{{kmer}}/tmp/"
    shell:
        """
        sleep 120
        rm -rf {params}
        echo 'cleaned up again'
        """
        