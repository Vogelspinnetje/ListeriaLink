import heapq
import math
import json


def top_five(seq: str, nor: int, avg_pps: float, header: str, subject: list[tuple[int, int, dict]], marker: int) -> list[tuple[int, int, dict]]:
    if len(subject) < 5:
        read_data = (marker, nor, {"read_length": len(seq),
                                         "phred": avg_pps,
                                         "header": header})
        heapq.heappush(subject, read_data)
    else:
        if marker > subject[0][0]:
            read_data = (marker, nor, {"read_length": len(seq),
                                             "phred": avg_pps,
                                             "header": header})
            heapq.heapreplace(subject, read_data)

    return subject


def distrubutions(marker: int, subject: dict) -> dict:
    afronding = str(marker)
    if afronding in subject:
        subject[afronding] += 1
    else:
        subject[afronding] = 1

    return subject


def read_fastq(fq_file_path: str, js_file_path: str) -> None:
    """
    Fungeert als de main voor de quality control
    """
    
    base_counts: dict = {"a":0,"c":0,"t":0,"g":0}
    number_of_reads: int = 0
    total_score: int = 0
    read_length_distribution: dict = {}
    phred_distribution: dict = {}
    longest_five_: list = []  
    best_five_: list = [] 

    with open(fq_file_path, 'r') as file:
        while True:
            header: str = file.readline().strip()
            if not header:
                break
            sequence: str= file.readline().strip()
            file.readline()
            quality:str = file.readline().strip()

            for basen in sequence:
                base_counts[basen.lower()] +=1

            number_of_reads +=1

            phred_scores: list[int] = [ord(char) - 33 for char in quality]
            p_values: list[float] = [10 ** (-score / 10) for score in phred_scores]
            total_score += sum(p_values)
            avg_phred_per_sequence: float = -10 * math.log10(sum(p_values) / len(
                sequence))

            longest_five_: list[tuple[int, int, dict]] = top_five(sequence, number_of_reads,
                                    avg_phred_per_sequence, header,
                                    longest_five_, len(sequence))

            best_five_: list[tuple[int, int, dict]] = top_five(sequence, number_of_reads,
                                    avg_phred_per_sequence, header,
                                    best_five_, avg_phred_per_sequence)

            read_length_distribution: dict = distrubutions(math.ceil((len(
                sequence)+1)/1000), read_length_distribution)

            phred_distribution: dict = distrubutions(int(round(
                avg_phred_per_sequence)), phred_distribution)


    total_bases: int = sum(base_counts.values())
    gc_ratio: float = (base_counts["g"]+base_counts["c"])/total_bases
    avg_phred: float = -10 * math.log10(total_score / total_bases)
    avg_length: float = total_bases / number_of_reads
    longest_five_: list[tuple[int, int, dict]] = sorted(longest_five_, key=lambda x: x[2]["read_length"],
                           reverse=True)
    best_five_: list[tuple[int, int, dict]] = sorted(best_five_, key=lambda x: x[2]["phred"], reverse=True)

    longest_five: list = []
    best_five: list = []

    for reads in longest_five_:
        longest_five.append(reads[2])
    for reads in best_five_:
        best_five.append(reads[2])

    data: dict = {"base_counts": base_counts,
            "gc_ratio": gc_ratio,
            "number_of_reads": number_of_reads,
            "total_bases": total_bases,
            "avg_phred": avg_phred,
            "avg_length": avg_length,
            "longest_five": longest_five,
            "best_five": best_five,
            "read_length_distribution": read_length_distribution,
            "phred_distribution": phred_distribution}

    with open(js_file_path, 'w') as file:
        json.dump(data, file, indent=4)


if __name__ == "__main__":
    fastq_file_path: str = snakemake.input.fastq
    json_file_path: str = snakemake.output.json

    read_fastq(fastq_file_path, json_file_path)