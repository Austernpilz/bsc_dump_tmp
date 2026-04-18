import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_kmer_counts(run_name, markov_txt, base_txt, random_txt, k_val, top_n=265):
    # Load and label
    df_m = pd.read_csv(markov_txt, sep='\t', names=['kmer', 'count']).assign(Type='Markov')
    df_b = pd.read_csv(base_txt, sep='\t', names=['kmer', 'count']).assign(Type='Reference')
    df_r = pd.read_csv(random_txt, sep='\t', names=['kmer', 'count']).assign(Type='Random')

    df_all = pd.concat([df_m, df_b, df_r]).sort_values('kmer')

    # Slice for readability if k is high
    kmers_to_show = df_all['kmer'].unique()[:top_n]
    df_plot = df_all[df_all['kmer'].isin(kmers_to_show)]

    plt.figure(figsize=(15, 7))
    sns.barplot(data=df_plot, x='kmer', y='count', hue='Type')
    plt.xticks(rotation=90, fontsize=9)
    plt.title(f'K-mer Counts for k={k_val} (Subset: {kmers_to_show[0]} to {kmers_to_show[-1]})')
    plt.ylabel('Frequency in Genome')
    plt.xlabel('K-mer sequence')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'results/plots/{run_name}_{k_val}kmer_plot.png')
    plt.close()


def plot_benchmarks(run_name, markov_bench, base_bench, random_bench, k_val, i="bwa|kmc"):
    """
    Compares Runtime (s) and Max Memory (MB) for the three mapping tasks.
    """
    paths = {'Markov': markov_bench, 'Reference': base_bench, 'Random': random_bench}
    data = []

    for label, path in paths.items():
        df = pd.read_csv(path, sep='\t')
        data.append({
            'Genome': label, 
            'Seconds': df.iloc[0]['s'], 
            'Memory_MB': df.iloc[0]['max_rss']
        })

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
    plt.savefig(f'results/plots/{run_name}_{k_val}_benchmark_plot_{i}.png')
    plt.close()

# RUN_NAME = config["name"]

rule plotting:
    input:
        txt=f"results/{RUN_NAME}/{{kmer}}/markov_kmc.txt",
        txt1=f"results/{RUN_NAME}/{{kmer}}/base_kmc.txt",
        txt2=f"results/{RUN_NAME}/{{kmer}}/random_kmc.txt",

        bwa=f"results/{RUN_NAME}/{{kmer}}/bwa-mem2.tsv",
        bwa1=f"results/{RUN_NAME}/base/bwa-mem2.tsv",
        bwa2=f"results/{RUN_NAME}/random/bwa-mem2.tsv",

        kmc=f"results/{RUN_NAME}/{{kmer}}/markov_kmc.tsv",
        kmc1=f"results/{RUN_NAME}/{{kmer}}/base_kmc.tsv",
        kmc2=f"results/{RUN_NAME}/{{kmer}}/random_kmc.tsv",

        msb=f"results/{RUN_NAME}/{{kmer}}/markov_sim.tsv",
        rsb=f"results/{RUN_NAME}/random/random_sim.tsv",


    output:
        f"results/plots/{RUN_NAME}_{{kmer}}kmer_plot.png",
        f"results/plots/{RUN_NAME}_{{kmer}}_benchmark_plot_bwa.png",
        f"results/plots/{RUN_NAME}_{{kmer}}_benchmark_plot_kmc.png",

    run:
        plot_kmer_counts(RUN_NAME, input.txt, input.txt1, input.txt2, wildcards.kmer)
        plot_benchmarks(RUN_NAME, input.bwa, input.bwa1, input.bwa2, wildcards.kmer, "bwa")
        plot_benchmarks(RUN_NAME, input.kmc, input.kmc1, input.kmc2, wildcards.kmer, "kmc")
        plot_benchmarks(RUN_NAME, input.msb, input.rsb, input.rsb, wildcards.kmer, "markov")

