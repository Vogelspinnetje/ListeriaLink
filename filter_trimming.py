import yaml
import math
from snakemake.script import snakemake


def read_fastq(file_path: str, fqe_path: str, tb: int, mp: int, ml: int) -> None:
    """
    Leest het file_path in. Controleert of alle parameters kloppen. Zo ja, dan wordt het weggeschreven naar fqe_path.
    Voer een FASTQ bestand in. FASTA werkt niet.
    """

    with open(file_path, 'r') as file, open(fqe_path, 'w') as outfile:
        while True:
            header: str = file.readline().strip()
            if not header:
                break
            sequence: str = file.readline().strip()
            plus: str = file.readline().strip()
            quality: str = file.readline().strip()

            trimmed_sequence: str = sequence[tb:-tb]
            trimmed_quality: str = quality[tb:-tb]

            if len(trimmed_sequence) < ml:
                continue

            phred_scores: list[int] = [ord(char) - 33 for char in quality]
            p_values: list[float] = [10 ** (-score / 10) for score in phred_scores]
            avg_phred_per_sequence: float = -10 * math.log10(sum(p_values) / len(
                sequence))

            if avg_phred_per_sequence >= mp:
                outfile.write(f"{header}\n{trimmed_sequence}\n"
                              f"{plus}\n{trimmed_quality}\n")


if __name__ == "__main__":
    with open('config.yaml') as yaml_file:
        config: dict = yaml.safe_load(yaml_file)
        
    TRIM_BASES: int = config['TRIM_BASES']
    MIN_PHRED: int = config['MIN_PHRED']
    MIN_LENGTH: int = config['MIN_LENGTH']
    fastq_file_path: str = snakemake.input.fastq
    fastq_edited_path: str = snakemake.output.fastq_edited

    read_fastq(fastq_file_path, fastq_edited_path, TRIM_BASES, MIN_PHRED,MIN_LENGTH)