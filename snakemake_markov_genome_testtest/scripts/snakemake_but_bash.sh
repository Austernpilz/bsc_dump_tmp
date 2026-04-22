


name: "covid"


#runs through different kmers for testing
kmer_start: 3
kmer_end: 27

COVID = "/bigdata/ag_abi/manuel/genome/covid/wuhCor1.fa"
C_bp = 29903
C_mb = 1
ECOLI = "/bigdata/ag_abi/manuel/genome/ecoli/GCF_003970875.1.fa"
E_bp = 4696000
E_mb = 5
FLY = "/bigdata/ag_abi/manuel/genome/fruit_fly/dm6.fa"
F_bp = 143726002
F_mb = 140
MOUS = "/bigdata/ag_abi/manuel/genome/mouse/mm39.fa"
M_bp = 2728222451
M_mb = 2600
HUMAN = "/bigdata/ag_abi/manuel/genome/human/hg38.fa"
H_bp = 3299210039
H_mb = 3200



run_markov_genome(){
    pixi run markov_genome \
        --input $1 \            #{input.fasta_file} 
        --output $2 \           #{output} 
        --order $3 \            #{params.k}
        --lens $4 #   {params.s} 2>> {log.err} 1>> {log.out} ;
        #scripts/fix_broken_files.sh {output} 
}

run_kmc(){
    kmc -fm -k$1 -t4  $2 $3 $4      #{input} {params.dump} {params.tmp}
    kmc_tools transform $3 dump $5
}

run_bwa(){
    bwa-mem2 mem -t16 $0 $1 > $2
    #{threads} {input.fasta_file} {input.reads} > {output.sam} \
}

run_random_markov(){
    markov_genome random 
}