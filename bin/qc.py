#!/usr/bin/env python3

import argparse
import csv
import shlex
import subprocess

import matplotlib.pyplot as plt
import pandas as pd

from Bio import SeqIO
from matplotlib.patches import Rectangle

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""


def make_qc_plot(
    depth_pos,
    n_density,
    ref_length,
    samplename,
    min_depth,
    amplicons,
    window=200,
    ylim_top=10**5,
    width_inches=20.0,
    height_inches=4.0,
):
    # Initialize position/depth to zero for every position.
    # Ensures that dataframe is correct dimensions. Otherwise script will fail for samples with zero depth across genome.
    depth_df = pd.DataFrame(
        {
            "position": [x + 1 for x in range(ref_length)],
            "depth": [0 for x in range(ref_length)],
        }
    )

    for pos in depth_pos:
        idx = int(pos[1]) - 1
        depth = int(pos[2])
        depth_df.loc[idx, "depth"] = depth

    depth_df["depth_moving_average"] = depth_df.iloc[:, 1].rolling(window=window).mean()

    n_df = pd.DataFrame(
        {
            "position": [pos[0] for pos in n_density],
            "n_density": [dens[1] for dens in n_density],
        }
    )

    fig, (ax_depth, ax_amplicons) = plt.subplots(
        2, gridspec_kw={"height_ratios": [5, 1]}
    )
    fig.set_size_inches(width_inches, height_inches)

    ax_n_density = ax_depth.twinx()

    ax_depth.set_xlabel("Position")

    ax_depth.set_ylabel("Depth", color="g")
    ax_depth.set_ylim(top=ylim_top, bottom=1)
    ax_depth.set_yscale("log")
    ax_depth.axhline(y=min_depth, c="blue", linestyle="dotted", linewidth=0.5)
    ax_depth.plot(depth_df["depth_moving_average"], color="g", linewidth=0.5)

    ax_n_density.set_ylabel("N density", color="r")
    ax_n_density.plot(n_df["n_density"], color="r", linewidth=0.5)
    ax_n_density.set_ylim(top=1, bottom=0)

    # ax_amplicons.get_shared_x_axes().join(ax_amplicons, ax_depth)
    ax_amplicons.sharex(ax_depth)
    ax_amplicons.plot()

    amplicon_pool_colors = {
        "1": "tab:red",
        "2": "tab:blue",
        "3": "tab:green",
        "4": "tab:orange",
    }

    for amplicon_num, amplicon in amplicons.items():
        amplicon_rectangle = Rectangle(
            (amplicon["start"], int(amplicon["pool"])),
            amplicon["length"],
            1,
            color=amplicon_pool_colors[amplicon["pool"]],
        )
        rx, ry = amplicon_rectangle.get_xy()
        cx = rx + amplicon_rectangle.get_width() / 2.0
        cy = ry + amplicon_rectangle.get_height() / 2.0
        ax_amplicons.add_patch(amplicon_rectangle)
        ax_amplicons.annotate(
            amplicon_num,
            (cx, cy),
            color="black",
            weight="bold",
            fontsize=6,
            ha="center",
            va="center",
        )

    ax_amplicons.set_ylim(top=4, bottom=0)
    ax_amplicons.yaxis.set_visible(False)
    ax_amplicons.xaxis.set_visible(False)

    plt.xlim(left=0)
    plt.title(samplename)
    plt.savefig(samplename + ".depth.png", bbox_inches="tight")


def read_depth_file(bamfile):
    p = subprocess.Popen(
        ["samtools", "depth", "-a", "-d", "0", bamfile], stdout=subprocess.PIPE
    )
    out, err = p.communicate()
    counter = 0

    pos_depth = []
    for ln in out.decode("utf-8").split("\n"):
        if ln:
            pos_depth.append(ln.split("\t"))

    return pos_depth


def read_primers(
    primer_bed_path,
    primer_name_delimiter="_",
    amplicon_number_offset=1,
):

    primers_by_name = {}
    with open(primer_bed_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            # I cannot believe that I have to do this, what sort of a god would allow such things to occur?
            # name = "_".join(str(x) for x in fields[3].split("_")[0:2])
            # pair_name = (
            #     name.replice("RIGHT", "LEFT")
            #     if "RIGHT" in name
            #     else name.replace("LEFT", "RIGHT")
            # )
            pool = fields[4]
            orientation = fields[5]
            amplicon_number = name.split(primer_name_delimiter)[amplicon_number_offset]
            primers_by_name[name] = {
                "name": name,
                "pair_name": "",
                "start": start,
                "end": end,
                "amplicon_number": amplicon_number,
                "pool": pool,
                "orientation": orientation,
            }

    for primer_name, primer in primers_by_name.items():
        for potential_pair_name, potential_pair_details in primers_by_name.items():
            if (
                potential_pair_details["amplicon_number"] == primer["amplicon_number"]
                and primer["orientation"] != potential_pair_details["orientation"]
            ):
                primers_by_name[primer_name]["pair_name"] = potential_pair_name
                primers_by_name[potential_pair_name]["pair_name"] = primer_name

    return primers_by_name


def primers_to_amplicons(primers):
    amplicons = {}
    positive_orientation_primers = [
        primers[p] for p in primers if primers[p]["orientation"] == "+"
    ]
    for primer in positive_orientation_primers:
        amplicon_number = primer["amplicon_number"]
        amplicon_pool = primer["pool"]
        amplicon_start = primer["end"]
        amplicon_end = primers[primer["pair_name"]]["start"]
        amplicon_length = amplicon_end - amplicon_start
        amplicons[amplicon_number] = {
            "number": amplicon_number,
            "start": amplicon_start,
            "end": amplicon_end,
            "length": amplicon_length,
            "pool": amplicon_pool,
        }

    return amplicons


def get_covered_pos(pos_depth, min_depth):
    counter = 0
    for contig, pos, depth in pos_depth:
        if int(depth) >= min_depth:
            counter = counter + 1

    return counter


def get_N_positions(fasta):
    n_pos = [i for i, letter in enumerate(fasta.seq.lower()) if letter == "n"]

    return n_pos


def get_pct_N_bases(fasta):

    count_N = len(get_N_positions(fasta))

    pct_N_bases = count_N / len(fasta.seq) * 100

    return pct_N_bases


def get_largest_N_gap(fasta):
    n_pos = get_N_positions(fasta)

    n_pos = [0] + n_pos + [len(fasta.seq)]

    n_gaps = [j - i for i, j in zip(n_pos[:-1], n_pos[1:])]

    return sorted(n_gaps)[-1]


def get_ref_length(ref):
    record = SeqIO.read(ref, "fasta")
    return len(record.seq)


def sliding_window_N_density(sequence, window=10):

    sliding_window_n_density = []
    for i in range(0, len(sequence.seq), 1):
        window_mid = i + (window / 2)
        window_seq = sequence.seq[i : i + window]
        n_count = window_seq.lower().count("n")
        n_density = n_count / window

        sliding_window_n_density.append([window_mid, n_density])

    return sliding_window_n_density


def get_num_reads(bamfile):

    st_filter = "0x900"
    command = "samtools view -c -F{} {}".format(st_filter, bamfile)
    what = shlex.split(command)

    return subprocess.check_output(what).decode().strip()


def main(args):
    primers = read_primers(args.primer_bed)

    with open(args.primer_pairs, "w") as f:
        for primer, data in primers.items():
            if data["orientation"] == "+":
                f.write(f"{data['name']}\t{data['pair_name']}\n")

    amplicons = primers_to_amplicons(primers)

    ## Depth calcs
    ref_length = get_ref_length(args.ref)
    depth_pos = read_depth_file(args.bam)

    depth_covered_bases = get_covered_pos(depth_pos, args.min_depth)

    pct_covered_bases = depth_covered_bases / ref_length * 100

    ## Number of aligned reads calculaton
    num_reads = get_num_reads(args.bam)

    # Unknown base calcs
    fasta = SeqIO.read(args.fasta, "fasta")

    pct_N_bases = 0
    largest_N_gap = 0
    qc_pass = "FALSE"

    if len(fasta.seq) != 0:

        pct_N_bases = get_pct_N_bases(fasta)
        largest_N_gap = get_largest_N_gap(fasta)

    qc_line = {
        "sample_name": args.sample,
        "pct_N_bases": "{:.2f}".format(pct_N_bases),
        "pct_covered_bases": "{:.2f}".format(pct_covered_bases),
        "longest_no_N_run": largest_N_gap,
        "num_aligned_reads": num_reads,
        "fasta": args.fasta,
        "bam": args.bam,
    }

    with open(args.outfile, "w") as csvfile:
        header = qc_line.keys()
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerow(qc_line)

    N_density = sliding_window_N_density(fasta)
    make_qc_plot(
        depth_pos,
        N_density,
        ref_length,
        args.sample,
        args.min_depth,
        amplicons,
        window=args.window,
        ylim_top=10**args.max_y_axis_exponent,
        width_inches=args.width,
        height_inches=args.height,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--primer-pairs", required=True)
    parser.add_argument("--primer-bed")
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--max-y-axis-exponent", type=int, default=4)
    parser.add_argument("--window", type=int, default=200)
    parser.add_argument("--width", type=float, default=36.0)
    parser.add_argument("--height", type=float, default=6.0)
    args = parser.parse_args()
    main(args)
