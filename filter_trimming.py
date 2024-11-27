import yaml
import math

from snakemake.script import snakemake


def read_fastq(file_path, fqe_path, tb, mp, ml):
    with open(file_path, 'r') as file, open(fqe_path, 'w') as outfile:
        while True:
            header = file.readline().strip()
            if not header:
                break
            sequence = file.readline().strip()
            plus = file.readline().strip()
            quality = file.readline().strip()

            trimmed_sequence = sequence[tb:-tb]
            trimmed_quality = quality[tb:-tb]

            if len(trimmed_sequence) < ml:
                continue

            phred_scores = [ord(char) - 33 for char in quality]
            p_values = [10 ** (-score / 10) for score in phred_scores]
            avg_phred_per_sequence = -10 * math.log10(sum(p_values) / len(
                sequence))

            if avg_phred_per_sequence >= mp:
                outfile.write(f"{header}\n{trimmed_sequence}\n"
                              f"{plus}\n{trimmed_quality}\n")


if __name__ == "__main__":
    with open('config.yaml') as yaml_file:
        config = yaml.safe_load(yaml_file)

    TRIM_BASES = config['TRIM_BASES']
    MIN_PHRED = config['MIN_PHRED']
    MIN_LENGTH = config['MIN_LENGTH']
    fastq_file_path = snakemake.input.fastq
    fastq_edited_path = snakemake.output.fastq_edited

    read_fastq(fastq_file_path, fastq_edited_path, TRIM_BASES, MIN_PHRED,MIN_LENGTH)