import os
import subprocess


def check_index(fasta_file, endings):
    if not os.path.isfile(fasta_file):
        print(f"{fasta_file} not a file")
        return False

    for suffix in endings:
        index_file = fasta_file+suffix

        if not os.path.isfile(index_file):
            print("{index_file} not a file")
            return False

        if os.path.getsize(index_file) == 0:
            print("size too small")
            return False

    return True

def clean_up(fasta_file, endings):
    for suffix in endings:
        index_file = fasta_file+suffix
        if os.path.isfile(index_file):
            os.remove(index_file)
            print(f"removed {index_file}")


def build_index(fasta_file):
    log_file = f"{fasta_file}_log.txt"
    try:
        with open(log_file, "w") as f:
            subprocess.run(
                f"bwa-mem2 index {fasta_file}",
                shell=True,
                check=True,
                stdout=f,
                stderr=f,
            )
        print(f"index build for {fasta_file}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error occurred! Exit code: {e.returncode}")
        print("___Error Details___")
        if os.path.isfile(log_file):
            with open(log_file, "r") as f:
                print("___Log Details___")
                print(f.read().splitlines()[-50:])
        return False

fasta_files = [
    "/bigdata/ag_abi/manuel/genome/covid/wuhCor1.fa",
    "/bigdata/ag_abi/manuel/genome/ecoli/GCF_003970875.1.fa",
    "/bigdata/ag_abi/manuel/genome/fruit_fly/dm6.fa",
    "/bigdata/ag_abi/manuel/genome/mouse/mm39.fa",
    "/bigdata/ag_abi/manuel/genome/human/hg38.fa",
]
endings = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]


for file in fasta_files:
    if not os.path.isfile(file):
        print(f"wrong file {file}")
    trys = 0
    while(trys<10):
        trys += 1
        if build_index(file) and check_index(file,endings):
            print(f"index build for {file}")
            break
        clean_up(file, endings)
