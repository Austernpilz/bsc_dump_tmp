import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import random
import subprocess

def kmer_to_int(kmer):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    # Calculate base-4 integer value
    val = 0
    for char in kmer:
        val = val * 4 + mapping.get(char, 0)
    return val

def plot_kmer_counts(output_name, names, markov_txt, base_txt, random_txt, k_val, top_n=1000):
    # Load and label
    df_m = pd.read_csv(markov_txt, sep='\t', names=['kmer', 'count']).assign(Type=names[0]).sort_values('kmer')
    df_b = pd.read_csv(base_txt, sep='\t', names=['kmer', 'count']).assign(Type=names[1]).sort_values('kmer')
    df_r = pd.read_csv(random_txt, sep='\t', names=['kmer', 'count']).assign(Type=names[2]).sort_values('kmer')

    df_all = pd.concat([df_m[:top_n], df_b[:top_n], df_r[:top_n]]).sort_values('kmer')

    kmers_to_show = df_all['kmer'].unique()[:top_n]

    df_plot = df_all[df_all['kmer'].isin(kmers_to_show)]
    df_plot.to_csv(output_name, sep='\t', index=False)

    df_plot['kmer_num'] = df_plot['kmer'].apply(kmer_to_int)

    plt.hist(df['kmer_num'], weights=df['count'], bins=100, alpha=0.5, label=df['Type'].iloc[0])

    plt.legend(loc='upper right')
    plt.title('Overlapping')

    plt.title(f'K-mer Counts for k={k_val}')
    plt.ylabel('Frequency in Genome')
    plt.xlabel('K-mer sequence')
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"results/{RUN_NAME}/{k_val}/output_name.png",)
    plt.close()



def plot_benchmarks(file_name, tsvv, k_val):
    """
    Compares Runtime (s) and Max Memory (MB) for the three mapping tasks.
    """
    data = []
    for p in tsvv:
        if os.path.isfile(p):
            data.append(p)

    for p in tsvv:
        if os.path.isfile(p):
            df = pd.read_csv(p, sep='\t')
            data.append({
                'Genome': os.path.basename(p), 
                'kmer' : k_val,
                'run' : RUN_NAME,
                'Seconds': df.iloc[0]['s'], 
                'Memory_MB': df.iloc[0]['max_rss']
            })
    if data:
        bench_df = pd.DataFrame(data)

        # Create a subplot with 1 row, 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Plot 1: Runtime
        sns.barplot(data=bench_df, x='Genome', y='Seconds', palette='magma', ax=ax1)
        ax1.set_title(f'Runtime (k={k_val})')
        ax1.set_ylabel('Seconds')

        # Plot 2: Memory
        sns.barplot(data=bench_df, x='Genome', y='Memory_MB', palette='viridis', ax=ax2)
        ax2.set_title(f'Peak Memory (k={k_val})')
        ax2.set_ylabel('Max RSS (MB)')

        plt.tight_layout()
        plt.savefig(file_name)
        plt.close()


rule plotting_kmc_kmer:
    input:
        f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
        markov = f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
        base = f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
        random = f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
        inter_markov = f"results/{RUN_NAME}/{{kmer}}/intersect_markov_base.txt",
        inter_random = f"results/{RUN_NAME}/{{kmer}}/intersect_random_base.txt",
    output:
        tsv = f"results/plots/kmer_{RUN_NAME}_{{kmer}}.tsv",
        inter0 = f"results/plots/kmer_{RUN_NAME}_{{kmer}}_intersect0.tsv",
        inter1 = f"results/plots/kmer_{RUN_NAME}_{{kmer}}_intersect1.tsv",
    run:
        plot_kmer_counts(output.tsv, ['random','base','markov'],input.markov, input.base, input.random, wildcards.kmer)
        plot_kmer_counts(output.inter0, ['markov_base','markov_random','base_random'], input.markov, input.base, input.inter_markov, wildcards.kmer)
        plot_kmer_counts(output.inter1, ['markov_base','markov_random','base_random'], input.random, input.base, input.inter_random, wildcards.kmer)


rule plotting_benchmark:
    input:
        f"results/{RUN_NAME}/random/random_sim.fa",
        f"results/{RUN_NAME}/random/sample.fastq",
        f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
        f"results/{RUN_NAME}/{{kmer}}/tmp/done.txt",
        f"results/{RUN_NAME}/{{qualifier}}/mapped.sam",
    output:
        f"results/plots/benchmark_markov_{RUN_NAME}_{{kmer}}kmer_plot.tsv",
    params:
        tsvvv = [ 
            f"results/{RUN_NAME}/base/bwa-mem2.tsv",
            f"results/{RUN_NAME}/random/bwa-mem2.tsv",
            f"results/{RUN_NAME}/{{kmer}}/bwa-mem2.tsv"
        ]
    run:
        plot_benchmarks(output, params.tsvvv, wildcards.kmer)

