import heapq
import math
import json


def top_five(seq, nor, avg_pps, header, subject, marker):
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


def distrubutions(marker, subject):
    afronding = str(marker)
    if afronding in subject:
        subject[afronding] += 1
    else:
        subject[afronding] = 1

    return subject


def read_fastq(fq_file_path, js_file_path):
    base_counts = {"a":0,"c":0,"t":0,"g":0}
    number_of_reads = 0
    total_score = 0
    read_length_distribution = {}
    phred_distribution = {}
    longest_five_ = []
    best_five_ = []

    with open(fq_file_path, 'r') as file:
        while True:
            header = file.readline().strip()
            if not header:
                break
            sequence = file.readline().strip()
            file.readline()
            quality = file.readline().strip()

            for basen in sequence:
                base_counts[basen.lower()] +=1

            number_of_reads +=1

            phred_scores = [ord(char) - 33 for char in quality]
            p_values = [10 ** (-score / 10) for score in phred_scores]
            total_score += sum(p_values)
            avg_phred_per_sequence = -10 * math.log10(sum(p_values) / len(
                sequence))

            longest_five_ = top_five(sequence, number_of_reads,
                                    avg_phred_per_sequence, header,
                                    longest_five_, len(sequence))

            best_five_ = top_five(sequence, number_of_reads,
                                    avg_phred_per_sequence, header,
                                    best_five_, avg_phred_per_sequence)

            read_length_distribution = distrubutions(math.ceil((len(
                sequence)+1)/1000), read_length_distribution)

            phred_distribution = distrubutions(int(round(
                avg_phred_per_sequence)), phred_distribution)


    total_bases = sum(base_counts.values())
    gc_ratio = (base_counts["g"]+base_counts["c"])/total_bases
    avg_phred = -10 * math.log10(total_score / total_bases)
    avg_length = total_bases / number_of_reads
    longest_five_ = sorted(longest_five_, key=lambda x: x[2]["read_length"],
                           reverse=True)
    best_five_ = sorted(best_five_, key=lambda x: x[2]["phred"], reverse=True)

    longest_five = []
    best_five = []
    for reads in longest_five_:
        longest_five.append(reads[2])
    for reads in best_five_:
        best_five.append(reads[2])

    data = {"base_counts": base_counts,
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

fastq_file_path = snakemake.input.fastq
json_file_path = snakemake.output.json
read_fastq(fastq_file_path, json_file_path)

