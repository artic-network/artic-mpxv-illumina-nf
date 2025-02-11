#!/usr/bin/env python

import os
import requests
import re
import json
import subprocess
import hashlib
import sys
import csv
import shutil

from clint.textui import colored
from Bio import SeqIO


def check_hash(filepath, manifest_hash):
    with open(filepath, "rb") as fh:
        data = fh.read()
        hash_md5 = hashlib.md5(data).hexdigest()
    if hash_md5 != manifest_hash:
        print(
            colored.yellow(
                f"md5 hash for {str(filepath)} does not match manifest, if this scheme has been downloaded previously, please ensure the file is not corrupted or has been tampered with. If this is a new download, please raise this as an issue here: https://github.com/quick-lab/primerschemes/issues"
            ),
            file=sys.stderr,
        )
        raise SystemExit(1)


def pick_best_ref(
    multi_ref_url: str,
    multi_ref_md5: str,
    read_file: str,
    n_reads: int,
    scheme_path: str,
    mm2_threads: int,
):
    """Pick the best reference from a multi-reference alignment

    Args:
        multi_ref_url (str): URL to the multi-reference alignment from quick-lab/primerschemes
        multi_ref_md5 (str): MD5 hash of the multi-reference alignment
        read_file (str): Path to the read file to test against the references
        n_reads (int): How many reads to sample from the read file for testing (default: 10000)
        scheme_path (str): Path to the scheme directory
        mm2_threads (int): Number of threads to use when aligning using minimap2

    Raises:
        SystemExit: If the ref selection fasta file cannot be fetched
        SystemExit: If the ref selection fasta file get returns a status code other than 200

    Returns:
        str: Primer scheme suffix of the most appropriate reference for the reads
    """

    ref_selection_path = os.path.join(scheme_path, "refselect.fasta")

    if not os.path.exists(scheme_path):
        os.makedirs(scheme_path, exist_ok=True)

    if not os.path.exists(ref_selection_path):
        try:
            response = requests.get(multi_ref_url)
        except requests.exceptions.RequestException as error:
            print(colored.red(f"Failed to ref selection fasta with Exception: {error}"))
            raise SystemExit(1)

        if response.status_code != 200:
            print(
                colored.red(
                    f"Failed to fetch ref selection fasta with status code: {response.status_code}"
                )
            )
            raise SystemExit(1)

        with open(ref_selection_path, "w") as ref_select_fasta:
            ref_select_fasta.write(response.text)

    check_hash(ref_selection_path, multi_ref_md5)

    flat_ref_path = os.path.join(scheme_path, "flattened_references.fasta")

    possible_references = {}

    with open(ref_selection_path, "r") as ref_selection_fh:

        for reference in SeqIO.parse(ref_selection_fh, "fasta"):
            suffix = reference.description.split()[1]
            possible_references[reference.id] = suffix

            # Flatten out the alignment into flat fasta reference
            flattened = str(reference.seq).replace("-", "")

            with open(flat_ref_path, "a") as flat_ref_fh:
                flat_ref_fh.write(f">{reference.description}\n{flattened}\n")

    # fq_it = mappy.fastx_read(fastq_path)
    seqtk = subprocess.run(
        ["seqtk", "sample", str(read_file), str(n_reads)], stdout=subprocess.PIPE
    )

    result = subprocess.run(
        [
            "minimap2",
            "-x",
            "map-ont",
            "-t",
            str(mm2_threads),
            str(flat_ref_path),
            "-",
        ],
        input=seqtk.stdout,
        stdout=subprocess.PIPE,
    )

    reader = csv.DictReader(
        result.stdout.decode("utf-8").split("\n"),
        delimiter="\t",
        fieldnames=[
            "query_name",
            "query_len",
            "query_start",
            "query_end",
            "strand",
            "ctg_name",
            "ctg_len",
            "ref_start",
            "ref_end",
            "n_matches",
            "alignment_len",
            "mapq",
        ],
    )

    read_results = {}

    for alignment in reader:

        if not alignment:
            continue

        identity = int(alignment["n_matches"]) / int(alignment["alignment_len"])

        read_results.setdefault(alignment["query_name"], {})
        read_results[alignment["query_name"]].setdefault(alignment["ctg_name"], 0)

        if identity > read_results[alignment["query_name"]][alignment["ctg_name"]]:
            read_results[alignment["query_name"]][alignment["ctg_name"]] = identity

    ref_results = {ref: 0 for ref in possible_references.keys()}

    for read, details in read_results.items():
        best_ctg = max(details, key=details.get)
        ref_results[best_ctg] += 1

    best_ref = max(ref_results, key=ref_results.get)

    return possible_references[best_ref]


def get_scheme(
    scheme_name: str,
    scheme_version: str,
    scheme_directory: str,
    scheme_length: int = False,
    read_file: str = False,
):
    """Get the primer scheme and reference fasta file from the manifest

    Args:
        scheme_name (str): Name of the scheme
        scheme_version (str): Version of the scheme
        scheme_directory (str): Directory where schemes are stored
        scheme_length (int, optional): Length of the scheme. Defaults to False.

    Returns:
        str: Path to the primer bed file
        str: Path to the reference fasta file
        str: Version of the scheme
    """

    try:
        response = requests.get(
            "https://raw.githubusercontent.com/quick-lab/primerschemes/main/index.json"
        )
        manifest = response.json()
        os.makedirs(scheme_directory, exist_ok=True)
        with open(os.path.join(scheme_directory, "index.json"), "w") as manifest_file:
            json.dump(manifest, manifest_file)

    except Exception:
        os.makedirs(scheme_directory, exist_ok=True)
        index_path = os.path.join(scheme_directory, "index.json")
        if os.path.exists(index_path):
            print(
                colored.yellow(
                    f"Could not get a remote copy of manifest but found a local one in {scheme_directory}, using this instead. Be warned this may be out of date"
                )
            )
            with open(index_path, "r") as manifest_file:
                manifest = json.load(manifest_file)
        else:
            print(
                colored.red(
                    f"Failed to fetch manifest and no local copy found in {scheme_directory}, please check your internet connection and try again"
                )
            )
            raise SystemExit(1)

    try:
        response = requests.get(
            "https://raw.githubusercontent.com/quick-lab/primerschemes/main/aliases.json"
        )
        aliases = response.json()
        os.makedirs(scheme_directory, exist_ok=True)
        with open(os.path.join(scheme_directory, "aliases.json"), "w") as aliases_file:
            json.dump(aliases, aliases_file)

    except Exception:
        os.makedirs(scheme_directory, exist_ok=True)
        aliases_path = os.path.join(scheme_directory, "aliases.json")
        if os.path.exists(aliases_path):
            print(
                colored.yellow(
                    f"Could not get a remote copy of manifest but found a local one in {scheme_directory}, using this instead. Be warned this may be out of date"
                )
            )
            with open(aliases_path, "r") as aliases_file:
                aliases = json.load(aliases_file)
        else:
            print(
                colored.red(
                    f"Failed to fetch manifest and no local copy found in {scheme_directory}, please check your internet connection and try again"
                )
            )
            raise SystemExit(1)

    scheme_name = scheme_name.lower().rstrip()

    if not manifest["primerschemes"].get(scheme_name):
        if not aliases.get(scheme_name):
            print(
                colored.red(
                    f"Requested scheme: {scheme_name} could not be found, please check https://github.com/quick-lab/primerschemes for available schemes or provide a scheme bed and fasta reference directly using --bed and --ref"
                )
            )
            raise SystemExit(1)

        scheme_name = aliases[scheme_name]

        if not manifest["primerschemes"].get(scheme_name):
            print(
                colored.red(
                    "Scheme name alias does not exist in manifest, this should never happen, please create an issue on https://github.com/quick-lab/primerschemes/issues or https://github.com/artic-network/fieldbioinformatics/issues if you see this message"
                )
            )
            raise SystemExit(1)

    scheme = manifest["primerschemes"][scheme_name]

    if len(scheme.keys()) > 1:
        if not scheme_length:
            print(
                colored.red(
                    f"Multiple lengths of the scheme {scheme_name} found, please provide a scheme length using --scheme-length, available lengths are: {', '.join(scheme.keys())}"
                )
            )
            raise SystemExit(1)

    else:
        scheme_length = list(scheme.keys())[0]

    if not scheme.get(scheme_length):
        print(
            colored.red(
                f"Provided length: {scheme_length} of the scheme {scheme_name} not found, please provide one of the following lengths using --scheme-length: {', '.join(scheme.keys())}"
            )
        )
        raise SystemExit(1)

    scheme_version = scheme_version.lower().rstrip()

    version_pattern = re.compile(r"v\d+\.\d+\.\d+")

    if not version_pattern.match(scheme_version):
        print(
            colored.red(
                "Invalid scheme version format, please provide a version in the format 'vX.X.X', e.g. v1.0.0"
            )
        )
        raise SystemExit(1)

    if not scheme[scheme_length].get(scheme_version):
        print(
            colored.red(
                f"Requested scheme version: {scheme_version} not found, available versions are: {', '.join(scheme[scheme_length].keys())}"
            )
        )
        raise SystemExit(1)

    if scheme[scheme_length][scheme_version].get("refselect"):
        if not read_file:
            print(
                colored.red(
                    f"Reference selection is available for this scheme but reads were not provided. Either provide a read file using --read-file or specify the reference to use by providing the same scheme name with the appropriate suffix, choices are: {', '.join(str(x) for x in scheme[scheme_length].keys() if '-' in x)}"
                )
            )
            raise SystemExit(1)

        print(
            colored.yellow(
                f"Reference selection is available for scheme {scheme_name}, deciding which reference to use based on your reads. If you would prefer to specify the reference to use, provide the same scheme name with the appropriate suffix, choices are: {', '.join(str(x) for x in scheme[scheme_length].keys() if '-' in x)}"
            )
        )

        if len(scheme[scheme_length][scheme_version]["refselect"].keys()) > 1:
            print(
                colored.red(
                    f"Multiple reference selection options found for scheme {scheme_name}, fieldbioinformatics does not currently support multi pathogen schemes"
                )
            )
            raise SystemExit(1)

        msa_data = next(
            iter(scheme[scheme_length][scheme_version]["refselect"].values())
        )

        scheme_select_dir = os.path.join(
            scheme_directory, scheme_name, scheme_length, scheme_version
        )

        suffix = pick_best_ref(
            multi_ref_url=msa_data["url"],
            multi_ref_md5=msa_data["md5"],
            read_file=read_file,
            n_reads=10000,
            scheme_path=scheme_select_dir,
            mm2_threads=4,
        )

        scheme_version = f"{scheme_version}-{suffix}"

        print(colored.yellow(f"Selected reference suffix: {suffix}"))

    scheme = scheme[scheme_length][scheme_version]

    scheme_path = os.path.join(
        scheme_directory, scheme_name, scheme_length, scheme_version
    )

    if not os.path.exists(scheme_path):
        os.makedirs(scheme_path, exist_ok=True)

    bed_url = scheme["primer_bed_url"]

    bed_name = bed_url.split("/")[-1]

    bed_path = os.path.join(scheme_path, bed_name)

    if not os.path.exists(bed_path):
        try:
            response = requests.get(bed_url)
        except requests.exceptions.RequestException as error:
            print(
                colored.red(
                    f"Failed to fetch primer bed file due to Exception: {error}"
                )
            )
            raise SystemExit(1)

        with open(bed_path, "w") as bed_file:
            bed_file.write(response.text)

    check_hash(bed_path, scheme["primer_bed_md5"])

    ref_url = scheme["reference_fasta_url"]

    ref_name = ref_url.split("/")[-1]

    ref_path = os.path.join(scheme_path, ref_name)

    if not os.path.exists(ref_path):
        try:
            response = requests.get(ref_url)
        except requests.exceptions.RequestException as error:
            print(
                colored.red(
                    f"Failed to fetch reference fasta file with Exception: {error}"
                )
            )
            raise SystemExit(1)

        if response.status_code != 200:
            print(
                colored.red(
                    f"Failed to fetch reference fasta file with status code: {response.status_code}"
                )
            )
            raise SystemExit(1)

        with open(ref_path, "w") as ref_file:
            ref_file.write(response.text)

    check_hash(ref_path, scheme["reference_fasta_md5"])

    return bed_path, ref_path, scheme_version


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Get a primer scheme and reference fasta file from the quick-lab/primerschemes repository"
    )
    parser.add_argument("--read-file", help="Read file to use for reference selection")
    parser.add_argument(
        "--scheme-directory",
        default="primerschemes",
        help="Directory to store primer schemes",
    )
    parser.add_argument(
        "scheme_name",
        help="Full name of scheme to get (as per the quick lab primerschemes repository e.g. 'artic-inrb-mpox/2500/v1.0.0')",
    )
    args = parser.parse_args()

    try:
        scheme_name, scheme_length, scheme_version = args.scheme_name.split("/")
    except ValueError:
        print(
            colored.red(
                f"Scheme name: {args.scheme_name} does not appear to conform to correct layout, please provide the full scheme name as per the quick lab primerschemes repository e.g. 'artic-inrb-mpox/2500/v1.0.0'"
            )
        )
        raise SystemExit(1)

    bed_path, ref_path, scheme_version = get_scheme(
        scheme_name=scheme_name,
        scheme_version=scheme_version,
        scheme_directory=args.scheme_directory,
        scheme_length=scheme_length,
        read_file=args.read_file,
    )

    os.makedirs("primer_scheme", exist_ok=True)

    shutil.copy(bed_path, "primer_scheme")
    shutil.copy(ref_path, "primer_scheme")


if __name__ == "__main__":
    main()
