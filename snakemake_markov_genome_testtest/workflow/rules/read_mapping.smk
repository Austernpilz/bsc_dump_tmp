rule bwa_mem2:
    input:
        fasta_file=config["ref_genom"],
        reads=f"results/{RUN_NAME}/{{qualifier}}/sample.fastq",
    output:
        sam = temp(f"results/{RUN_NAME}/{{qualifier}}/mapped.sam"),
        done = f"results/{RUN_NAME}/{{qualifier}}/mapped.txt",
    log:
        out=f"results/{RUN_NAME}/{{qualifier}}/bwa-mem2.log"
        err=f"results/{RUN_NAME}/{{qualifier}}/bwa-mem2.log"
    benchmark:
        f"results/{RUN_NAME}/{{qualifier}}/bwa-mem2.tsv"
    wildcard_constraints:
        qualifier="base|random|[0-9]+"
    priority: 
        100
    threads: config["bwa_t"]
    resources: mem_mb=s=config["mb"]* (28 + config["bwa_t"])
    shell:
        """
        bwa-mem2 mem -t {threads} {input.fasta_file} {input.reads} > {output.sam} \
        2>> {log.err} 1>> {log.out}
        touch {output.done}
        """