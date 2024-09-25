#!/usr/bin/env python

from copy import copy
import csv
import sys
import re
import random
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]

# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]

def merge_sites(canonical, alt):
    """Merges a canonical primer site with an alt site, producing an interval that encompasses both

    Parameters
    ----------
    canonical : dict
        The canonical primer site, provided as a dictionary of the bed file row
    alt : dict
        The alt primer site, provided as a dictionary of the bed file row

    Returns
    -------
    dict
        A dictionary of the merged site, where the dict represents a bed file row
    """
    # base the merged site on the canonical
    mergedSite = canonical

    # check the both the canonical and alt are the same direction
    if canonical["direction"] != alt["direction"]:
        print(
            "could not merge alt with different orientation to canonical",
            file=sys.stderr,
        )
        raise SystemExit(1)

    # merge the start/ends of the alt with the canonical to get the largest window possible
    if alt["start"] < canonical["start"]:
        mergedSite["start"] = alt["start"]
    if alt["end"] > canonical["end"]:
        mergedSite["end"] = alt["end"]
    return mergedSite


def identify_bed_file(bed_file):
    """Identifies the version of the primer ID format used in the bed file

    Parameters
    ----------
    bed_file : str
        The bed file to identify the primer ID format

    Returns
    -------
    int
        The version of the primer ID format used in the bed file
    """

    V1_pattern = re.compile(
        r"[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)(_ALT[0-9]*|_alt[0-9]*)*"
    )

    V2_pattern = re.compile(r"[a-zA-Z0-9\-]+_[0-9]+_(LEFT|RIGHT)_[0-9]+")

    version = False

    with open(bed_file, "r") as bed_fh:
        bed = bed_fh.readlines()

    for line in bed:
        if line.startswith("#"):
            continue

        splits = line.strip().split("\t")

        if len(splits) < 6:
            print(
                    f"Invalid bed file format, only found {len(splits)} columns. For valid formats please see https://github.com/quick-lab/primerschemes"
            )
            raise SystemExit(1)

        if not V1_pattern.match(splits[3]) and not V2_pattern.match(splits[3]):
            print(
                    f"Invalid primer ID format, {splits[3]} does not match the expected format"
            )
            raise SystemExit(1)

        if V2_pattern.match(splits[3]):
            if len(splits) < 7:
                print(
                        f"Invalid bed file format, only found {len(splits)} columns. For valid formats please see https://github.com/ChrisgKent/primal-page"
                )

            if not version:
                version = 3

            if version != 3:
                print(
                        f"Scheme BED does not appear to be a consistent scheme version, for primer scheme formats please see https://github.com/ChrisgKent/primal-page"
                )
                raise SystemExit(1)

        elif V1_pattern.match(splits[3]):

            if len(splits) == 7:
                if not version:
                    version = 2
            elif len(splits) == 6:
                if not version:
                    version = 1
            else:
                print(
                        f"Invalid bed file format, found {len(splits)} columns with V1 primer names. For valid formats please see https://github.com/ChrisgKent/primal-page"
                )

            if version == 3:
                print(
                        f"Scheme BED mixed primer ID formats, please ensure your BED file is consistent"
                )
                raise SystemExit(1)

        return version


def read_bed_file(fn):
    """Parses a V2/V3 bed file and collapses primers into canonical primer sites

    Parameters
    ----------
    fn : str
        The bedfile to parse

    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    version = identify_bed_file(fn)

    if version in (1, 2):
        return read_bed_file_legacy(fn)

    primers = pd.read_csv(
        fn,
        sep="\t",
        comment="#",
        header=None,
        names=["chrom", "start", "end", "Primer_ID", "PoolName", "direction"],
        dtype={
            "chrom": str,
            "start": int,
            "end": int,
            "Primer_ID": str,
            "PoolName": str,
        },
        usecols=(0, 1, 2, 3, 4, 5),
        skiprows=0,
    )
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    canonical_primers = {}
    for _, row in primers.iterrows():
        scheme_name, primer_id, direction, primer_n = row["Primer_ID"].split("_")

        if (primer_id, direction) not in canonical_primers:
            canonical_primers[(primer_id, direction)] = row.to_dict()
            continue

        canonical_primers[(primer_id, direction)] = merge_sites(
            canonical_primers[(primer_id, direction)], row
        )

    # return the bedFile as a list
    return [value for value in canonical_primers.values()]


def read_bed_file_legacy(fn):
    """Parses a bed file and collapses alts into canonical primer sites

    Parameters
    ----------
    fn : str
        The bedfile to parse

    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    primers = pd.read_csv(
        fn,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "Primer_ID", "PoolName", "direction"],
        dtype={
            "chrom": str,
            "start": int,
            "end": int,
            "Primer_ID": str,
            "PoolName": str,
        },
        usecols=(0, 1, 2, 3, 4, 5),
        skiprows=0,
    )
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    # separate alt primers into a new dataframe
    altFilter = primers["Primer_ID"].str.contains("_alt")
    alts = pd.DataFrame(
        columns=("chrom", "start", "end", "Primer_ID", "PoolName", "direction")
    )
    alts = pd.concat([alts, primers[altFilter]])
    primers = primers.drop(primers[altFilter].index.values)

    # convert the primers dataframe to dictionary, indexed by Primer_ID
    #  - verify_integrity is used to prevent duplicate Primer_IDs being processed
    bedFile = primers.set_index(
        "Primer_ID", drop=False, verify_integrity=True
    ).T.to_dict()

    # if there were no alts, return the bedfile as a list of dicts
    if len(alts.index) == 0:
        return list(bedFile.values())

    # merge alts
    for _, row in alts.iterrows():
        primerID = row["Primer_ID"].split("_alt")[0]

        # check the bedFile if another version of this primer exists
        if primerID not in bedFile:

            # add to the bed file and continue
            bedFile[primerID] = row.to_dict()
            continue

        # otherwise, we've got a primer ID we've already seen so merge the alt
        mergedSite = merge_sites(bedFile[primerID], row)

        # update the bedFile
        bedFile[primerID] = mergedSite

    # return the bedFile as a list
    return [value for value in bedFile.values()]

def find_primer(bed, pos, direction, threshold=20):
    """Given a reference position and a direction of travel, walk out and find the nearest primer site.

    Parameters
    ----------
    bed : list
        A list of dictionaries, where each dictionary contains a row of bedfile data
    pos : int
        The position in the reference sequence to start from
    direction : string
        The direction to search along the reference sequence

    Returns
    -------
    tuple[int, int, dict] | bool
        A tuple containing the distance to the primer, the relative position of the primer, and the primer site, or False if no primer found
    """
    from operator import itemgetter

    if direction == "+":
        primer_distances = [
            (abs(p["start"] - pos), p["start"] - pos, p)
            for p in bed
            if (p["direction"] == direction) and (pos >= (p["start"] - threshold))
        ]

    else:
        primer_distances = [
            (abs(p["end"] - pos), p["end"] - pos, p)
            for p in bed
            if (p["direction"] == direction) and (pos <= (p["end"] + threshold))
        ]

    if not primer_distances:
        return False

    closest = min(
        primer_distances,
        key=itemgetter(0),
    )

    return closest


def trim(segment, primer_pos, end, debug):
    """Soft mask an alignment to fit within primer start/end sites.

    Parameters
    ----------
    segment : pysam.AlignedSegment
        The aligned segment to mask
    primer_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    debug : bool
        If True, will print soft masking info during trimming
    """
    # get a copy of the cigar tuples to work with
    cigar = copy(segment.cigartuples)

    # get the segment position in the reference (depends on if start or end of the segment is being processed)
    if not end:
        pos = segment.pos
    else:
        pos = segment.reference_end

    # process the CIGAR to determine how much softmasking is required
    eaten = 0
    while 1:

        # chomp CIGAR operations from the start/end of the CIGAR
        try:
            if end:
                flag, length = cigar.pop()
            else:
                flag, length = cigar.pop(0)
            if debug:
                print("Chomped a %s, %s" % (flag, length), file=sys.stderr)
        except IndexError:
            print(
                "Ran out of cigar during soft masking - completely masked read will be ignored",
                file=sys.stderr,
            )
            break

        # if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if consumesReference[flag]:
            if not end:
                pos += length
            else:
                pos -= length

        # if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if consumesQuery[flag]:
            eaten += length

        # stop processing the CIGAR if we've gone far enough to mask the primer
        if not end and pos >= primer_pos and flag == 0:
            break
        if end and pos <= primer_pos and flag == 0:
            break

    # calculate how many extra matches are needed in the CIGAR
    extra = abs(pos - primer_pos)
    if debug:
        print("extra %s" % (extra), file=sys.stderr)
    if extra:
        if debug:
            print("Inserted a %s, %s" % (0, extra), file=sys.stderr)
        if end:
            cigar.append((0, extra))
        else:
            cigar.insert(0, (0, extra))
        eaten -= extra

    # softmask the left primer
    if not end:

        # update the position of the leftmost mappinng base
        segment.pos = pos - extra
        if debug:
            print("New pos: %s" % (segment.pos), file=sys.stderr)

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if debug:
                print(
                    "softmask created a leading deletion in the CIGAR, shuffling the alignment",
                    file=sys.stderr,
                )
            while 1:
                if cigar[0][0] != 2:
                    break
                _, length = cigar.pop(0)
                segment.pos += length

        # now add the leading softmask
        cigar.insert(0, (4, eaten))

    # softmask the right primer
    else:
        cigar.append((4, eaten))

    # check the new CIGAR and replace the old one
    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        raise ("invalid cigar operation created - possibly due to INDEL in primer")
    segment.cigartuples = cigar
    return segment

def handle_paired_segment(
    segments: tuple [pysam.AlignedSegment, pysam.AlignedSegment],
    bed: dict,
    args: argparse.Namespace,
    min_mapq: int,
    report_writer: csv.DictWriter = False,
):
    """Handle the alignment segment including

    Args:
        segment (pysam.AlignedSegment): The alignment segment to process
        bed (dict): The primer scheme
        reportfh (typing.IO): The report file handle
        args (argparse.Namespace): The command line arguments

    Returns:
        tuple [int, pysam.AlignedSegment] | bool: A tuple containing the amplicon number and the alignment segment, or False if the segment is to be skipped
    """

    segment1, segment2 = segments

    if not segment1 or not segment2:
        print("Segment pair skipped as at least one segment in pair does not exist", file=sys.stderr)
        return False

    # filter out unmapped and supplementary alignment segments
    if segment1.is_unmapped or segment2.is_unmapped:
        print("Segment pair: %s skipped as unmapped" % (segment1.query_name), file=sys.stderr)
        return False
    if segment1.is_supplementary or segment2.is_supplementary:
        print("Segment pair: %s skipped as supplementary" % (segment1.query_name), file=sys.stderr)
        return False
    if segment1.mapping_quality < min_mapq or segment2.mapping_quality < min_mapq:
        print(
            "Segment pair: %s skipped as mapping quality below threshold" % (segment1.query_name),
            file=sys.stderr,
        )
        return False

    # locate the nearest primers to this alignment segment
    p1 = find_primer(bed, segment1.reference_start, "+", args.primer_match_threshold)
    p2 = find_primer(bed, segment2.reference_end, "-", args.primer_match_threshold)

    if not p1 or not p2:
        print(
            "Paired segment: %s skipped as no primer found for segment" % (segment1.query_name),
            file=sys.stderr,
        )
        return False

    # check if primers are correctly paired and then assign read group
    # NOTE: removed this as a function as only called once
    # TODO: will try improving this / moving it to the primer scheme processing code
    correctly_paired = p1[2]["Primer_ID"].replace("_LEFT", "") == p2[2][
        "Primer_ID"
    ].replace("_RIGHT", "")

    if not args.no_read_groups:
        if correctly_paired:
            segment1.set_tag("RG", p1[2]["PoolName"])
            segment2.set_tag("RG", p1[2]["PoolName"])
        else:
            segment1.set_tag("RG", "unmatched")
            segment2.set_tag("RG", "unmatched")

    # get the amplicon number
    amplicon = p1[2]["Primer_ID"].split("_")[1]

    if args.report:
        # update the report with this alignment segment + primer details
        report = {
            "QueryName": segment1.query_name,
            "ReferenceStart": segment1.reference_start,
            "ReferenceEnd": segment2.reference_end,
            "PrimerPair": f"{p1[2]['Primer_ID']}_{p2[2]['Primer_ID']}",
            "Primer1": p1[2]["Primer_ID"],
            "Primer1Start": abs(p1[1]),
            "Primer2": p2[2]["Primer_ID"],
            "Primer2Start": abs(p2[1]),
            "IsSecondary": segment1.is_secondary,
            "IsSupplementary": segment1.is_supplementary,
            "Start": p1[2]["start"],
            "End": p2[2]["end"],
            "CorrectlyPaired": correctly_paired,
        }

        report_writer.writerow(report)

    if args.remove_incorrect_pairs and not correctly_paired:
        print(
            "Paired segment: %s skipped as not correctly paired" % (segment1.query_name),
            file=sys.stderr,
        )
        return False

    if args.verbose:
        # Dont screw with the order of the dict
        report_str = "\t".join(str(x) for x in report.values())
        print(report_str, file=sys.stderr)

    # get the primer positions
    if args.trim_primers:
        p1_position = p1[2]["end"]
        p2_position = p2[2]["start"]
    else:
        p1_position = p1[2]["start"]
        p2_position = p2[2]["end"]

    # softmask the alignment if left primer start/end inside alignment
    if segment1.reference_start < p1_position:
        try:
            trim(segment1, p1_position, False, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment1.reference_start, p1_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking left primer in {} (error: {}), skipping".format(
                    segment1.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # softmask the alignment if right primer start/end inside alignment
    if segment2.reference_end > p2_position:
        try:
            trim(segment2, p2_position, True, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment2.reference_start, p2_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking right primer in {} (error: {}), skipping".format(
                    segment1.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # check the the alignment still contains bases matching the reference
    if "M" not in segment1.cigarstring or "M" not in segment2.cigarstring:
        print(
            "Paired segment: %s dropped as does not match reference post masking"
            % (segment1.query_name),
            file=sys.stderr,
        )
        return False

    return (amplicon, segments)

def handle_segment(
    segment: pysam.AlignedSegment,
    bed: dict,
    args: argparse.Namespace,
    min_mapq: int,
    report_writer: csv.DictWriter = False,
):
    """Handle the alignment segment including

    Args:
        segment (pysam.AlignedSegment): The alignment segment to process
        bed (dict): The primer scheme
        reportfh (typing.IO): The report file handle
        args (argparse.Namespace): The command line arguments

    Returns:
        tuple [int, pysam.AlignedSegment] | bool: A tuple containing the amplicon number and the alignment segment, or False if the segment is to be skipped
    """

    # filter out unmapped and supplementary alignment segments
    if segment.is_unmapped:
        print("%s skipped as unmapped" % (segment.query_name), file=sys.stderr)
        return False
    if segment.is_supplementary:
        print("%s skipped as supplementary" % (segment.query_name), file=sys.stderr)
        return False
    if segment.mapping_quality < min_mapq:
        print(
            "%s skipped as mapping quality below threshold" % (segment.query_name),
            file=sys.stderr,
        )
        return False

    # locate the nearest primers to this alignment segment
    p1 = find_primer(bed, segment.reference_start, "+", args.primer_match_threshold)
    p2 = find_primer(bed, segment.reference_end, "-", args.primer_match_threshold)

    if not p1 or not p2:
        print(
            "%s skipped as no primer found for segment" % (segment.query_name),
            file=sys.stderr,
        )
        return False

    # check if primers are correctly paired and then assign read group
    # NOTE: removed this as a function as only called once
    # TODO: will try improving this / moving it to the primer scheme processing code
    correctly_paired = p1[2]["Primer_ID"].replace("_LEFT", "") == p2[2][
        "Primer_ID"
    ].replace("_RIGHT", "")

    if not args.no_read_groups:
        if correctly_paired:
            segment.set_tag("RG", p1[2]["PoolName"])
        else:
            segment.set_tag("RG", "unmatched")

    # get the amplicon number
    amplicon = p1[2]["Primer_ID"].split("_")[1]

    if args.report:
        # update the report with this alignment segment + primer details
        report = {
            "QueryName": segment.query_name,
            "ReferenceStart": segment.reference_start,
            "ReferenceEnd": segment.reference_end,
            "PrimerPair": f"{p1[2]['Primer_ID']}_{p2[2]['Primer_ID']}",
            "Primer1": p1[2]["Primer_ID"],
            "Primer1Start": abs(p1[1]),
            "Primer2": p2[2]["Primer_ID"],
            "Primer2Start": abs(p2[1]),
            "IsSecondary": segment.is_secondary,
            "IsSupplementary": segment.is_supplementary,
            "Start": p1[2]["start"],
            "End": p2[2]["end"],
            "CorrectlyPaired": correctly_paired,
        }

        report_writer.writerow(report)

    if args.remove_incorrect_pairs and not correctly_paired:
        print(
            "%s skipped as not correctly paired" % (segment.query_name),
            file=sys.stderr,
        )
        return False

    if args.verbose:
        # Dont screw with the order of the dict
        report_str = "\t".join(str(x) for x in report.values())
        print(report_str, file=sys.stderr)

    # get the primer positions
    if args.trim_primers:
        p1_position = p1[2]["end"]
        p2_position = p2[2]["start"]
    else:
        p1_position = p1[2]["start"]
        p2_position = p2[2]["end"]

    # softmask the alignment if left primer start/end inside alignment
    if segment.reference_start < p1_position:
        try:
            trim(segment, p1_position, False, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment.reference_start, p1_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking left primer in {} (error: {}), skipping".format(
                    segment.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # softmask the alignment if right primer start/end inside alignment
    if segment.reference_end > p2_position:
        try:
            trim(segment, p2_position, True, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment.reference_start, p2_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking right primer in {} (error: {}), skipping".format(
                    segment.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # check the the alignment still contains bases matching the reference
    if "M" not in segment.cigarstring:
        print(
            "%s dropped as does not match reference post masking"
            % (segment.query_name),
            file=sys.stderr,
        )
        return False

    return (amplicon, segment)


def generate_amplicons(bed: list):
    """Generate a dictionary of amplicons from a primer scheme list (generated by vcftagprimersites/read_bed_file)

    Args:
        bed (list): A list of dictionaries, where each dictionary contains a row of bedfile data (generated by vcftagprimersites/read_bed_file), assumes that all redundant primers have been expanded

    Raises:
        ValueError: Primer direction not recognised

    Returns:
        dict: A dictionary of amplicons, where each key is the amplicon number and the value is a dictionary containing the primer start, primer end, insert start, insert end, length and circularity
    """

    amplicons = {}
    for primer in bed:

        amplicon = primer["Primer_ID"].split("_")[1]

        amplicons.setdefault(amplicon, {})

        if primer["direction"] == "+":
            amplicons[amplicon]["p_start"] = primer["start"]
            amplicons[amplicon]["start"] = primer["end"] + 1

        elif primer["direction"] == "-":
            amplicons[amplicon]["p_end"] = primer["end"]
            amplicons[amplicon]["end"] = primer["start"] - 1

        else:
            raise ValueError("Primer direction not recognised")

    for amplicon in amplicons:
        if not all([x in amplicons[amplicon] for x in ["p_start", "p_end"]]):
            raise ValueError(f"Primer scheme for amplicon {amplicon} is incomplete")

        # Check if primer runs accross reference start / end -> circular virus
        amplicons[amplicon]["circular"] = (
            amplicons[amplicon]["p_start"] > amplicons[amplicon]["p_end"]
        )

        # Calculate amplicon length considering that the "length" may be negative if the genome is circular
        amplicons[amplicon]["length"] = abs(
            amplicons[amplicon]["p_end"] - amplicons[amplicon]["p_start"]
        )

    return amplicons

def normalise_paired(trimmed_segments: dict, normalise: int, bed: list):
    """Normalise the depth of the trimmed segments to a given value. Perform per-amplicon normalisation using numpy vector maths to determine whether the segment in question would take the depth closer to the desired depth accross the amplicon.

    Args:
        trimmed_segments (dict): Dict containing amplicon number as key and list of tuples liek: [pysam.AlignedSegment, pysam.AlignedSegment] as value
        normalise (int): Desired normalised depth
        bed (list): Primer scheme list (generated by vcftagprimersites/read_bed_file)
        trim_primers (bool): Whether to trim primers from the reads

    Raises:
        ValueError: Amplicon assigned to segment not found in primer scheme file

    Returns:
       list : List of pysam.AlignedSegment to output
    """

    amplicons = generate_amplicons(bed)

    output_segments = []

    mean_depths = {x: 0 for x in amplicons}

    for amplicon, segments in trimmed_segments.items():
        if amplicon not in amplicons:
            raise ValueError(f"Segment {amplicon} not found in primer scheme file")

        desired_depth = np.full_like(
            (amplicons[amplicon]["length"],), normalise, dtype=int
        )

        amplicon_depth = np.zeros((amplicons[amplicon]["length"],), dtype=int)

        if not segments:
            print(
                f"No segments assigned to amplicon {amplicon}, skipping",
                file=sys.stderr,
            )
            continue

        random.shuffle(segments)

        distance = np.mean(np.abs(amplicon_depth - desired_depth))

        for paired_segments in segments:

            test_depths = np.copy(amplicon_depth)

            segment1, segment2 = paired_segments

            for segment in (segment1, segment2):

                relative_start = segment.reference_start - amplicons[amplicon]["p_start"]

                if relative_start < 0:
                    relative_start = 0

                relative_end = segment.reference_end - amplicons[amplicon]["p_start"]

                test_depths[relative_start:relative_end] += 1

            test_distance = np.mean(np.abs(test_depths - desired_depth))

            if test_distance < distance:
                amplicon_depth = test_depths
                distance = test_distance
                output_segments.append(segment1)
                output_segments.append(segment2)

        mean_depths[amplicon] = np.mean(amplicon_depth)

    return output_segments, mean_depths

def normalise(trimmed_segments: dict, normalise: int, bed: list):
    """Normalise the depth of the trimmed segments to a given value. Perform per-amplicon normalisation using numpy vector maths to determine whether the segment in question would take the depth closer to the desired depth accross the amplicon.

    Args:
        trimmed_segments (dict): Dict containing amplicon number as key and list of pysam.AlignedSegment as value
        normalise (int): Desired normalised depth
        bed (list): Primer scheme list (generated by vcftagprimersites/read_bed_file)
        trim_primers (bool): Whether to trim primers from the reads

    Raises:
        ValueError: Amplicon assigned to segment not found in primer scheme file

    Returns:
       list : List of pysam.AlignedSegment to output
    """

    amplicons = generate_amplicons(bed)

    output_segments = []

    mean_depths = {x: 0 for x in amplicons}

    for amplicon, segments in trimmed_segments.items():
        if amplicon not in amplicons:
            raise ValueError(f"Segment {amplicon} not found in primer scheme file")

        desired_depth = np.full_like(
            (amplicons[amplicon]["length"],), normalise, dtype=int
        )

        amplicon_depth = np.zeros((amplicons[amplicon]["length"],), dtype=int)

        if not segments:
            print(
                f"No segments assigned to amplicon {amplicon}, skipping",
                file=sys.stderr,
            )
            continue

        random.shuffle(segments)

        distance = np.mean(np.abs(amplicon_depth - desired_depth))

        for segment in segments:
            test_depths = np.copy(amplicon_depth)

            relative_start = segment.reference_start - amplicons[amplicon]["p_start"]

            if relative_start < 0:
                relative_start = 0

            relative_end = segment.reference_end - amplicons[amplicon]["p_start"]

            test_depths[relative_start:relative_end] += 1

            test_distance = np.mean(np.abs(test_depths - desired_depth))

            if test_distance < distance:
                amplicon_depth = test_depths
                distance = test_distance
                output_segments.append(segment)

        mean_depths[amplicon] = np.mean(amplicon_depth)

    return output_segments, mean_depths

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not read.is_proper_pair:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def go(args):
    """Filter and soft mask an alignment file so that the alignment boundaries match the primer start and end sites.

    Based on the most likely primer position, based on the alignment coordinates.
    """
    # prepare the report outfile
    if args.report:
        reportfh = open(args.report, "w")
        report_headers = [
            "QueryName",
            "ReferenceStart",
            "ReferenceEnd",
            "PrimerPair",
            "Primer1",
            "Primer1Start",
            "Primer2",
            "Primer2Start",
            "IsSecondary",
            "IsSupplementary",
            "Start",
            "End",
            "CorrectlyPaired",
        ]
        report_writer = csv.DictWriter(reportfh, fieldnames=report_headers)
        report_writer.writeheader()

    # open the primer scheme and get the pools
    bed = read_bed_file(args.bedfile)
    pools = set([row["PoolName"] for row in bed])
    pools.add("unmatched")

    # open the input SAM file and process read groups
    infile = pysam.AlignmentFile("-", "rb")
    bam_header = infile.header.copy().to_dict()
    if not args.no_read_groups:
        bam_header["RG"] = []
        for pool in pools:
            read_group = {}
            read_group["ID"] = pool
            bam_header["RG"].append(read_group)

    # prepare the alignment outfile
    outfile = pysam.AlignmentFile("-", "wh", header=bam_header)

    trimmed_segments = dict()

    if args.paired:
        read_pairs = read_pair_generator(infile)

        for segments in read_pairs:
            if args.report:
                trimming_tuple = handle_paired_segment(
                    segments=segments,
                    bed=bed,
                    args=args,
                    report_writer=report_writer,
                    min_mapq=args.min_mapq,
                )
            else:
                trimming_tuple = handle_paired_segment(
                    segments=segments,
                    bed=bed,
                    args=args,
                    min_mapq=args.min_mapq,
                )
            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_pair = trimming_tuple
            trimmed_segments.setdefault(amplicon, [])

            if trimmed_segments:
                trimmed_segments[amplicon].append(trimmed_pair)

        if args.normalise:
            output_segments, mean_amp_depths = normalise_paired(
                trimmed_segments, args.normalise, bed
            )

            # write mean amplicon depths to file
            if args.amp_depth_report:
                with open(args.amp_depth_report, "w") as amp_depth_report_fh:
                    amp_depth_report_fh.write("amplicon\tmean_depth\n")
                    for amplicon, depth in mean_amp_depths.items():
                        amp_depth_report_fh.write(f"{amplicon}\t{depth}\n")

            for output_segment in output_segments:
                outfile.write(output_segment)

    else:
        # iterate over the alignment segments in the input SAM file
        for segment in infile:
            if args.report:
                trimming_tuple = handle_segment(
                    segment=segment,
                    bed=bed,
                    args=args,
                    report_writer=report_writer,
                    min_mapq=args.min_mapq,
                )
            else:
                trimming_tuple = handle_segment(
                    segment=segment,
                    bed=bed,
                    args=args,
                    min_mapq=args.min_mapq,
                )
            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_segment = trimming_tuple
            trimmed_segments.setdefault(amplicon, [])

            if trimmed_segment:
                trimmed_segments[amplicon].append(trimmed_segment)

        # normalise if requested
        if args.normalise:
            output_segments, mean_amp_depths = normalise(
                trimmed_segments, args.normalise, bed
            )

            # write mean amplicon depths to file
            if args.amp_depth_report:
                with open(args.amp_depth_report, "w") as amp_depth_report_fh:
                    amp_depth_report_fh.write("amplicon\tmean_depth\n")
                    for amplicon, depth in mean_amp_depths.items():
                        amp_depth_report_fh.write(f"{amplicon}\t{depth}\n")

            for output_segment in output_segments:
                outfile.write(output_segment)
        else:
            for amplicon, segments in trimmed_segments.items():
                for segment in segments:
                    outfile.write(segment)

    # close up the file handles
    infile.close()
    outfile.close()
    if args.report:
        reportfh.close()


def main():
    parser = argparse.ArgumentParser(
        description="Trim alignments from an amplicon scheme."
    )
    parser.add_argument("bedfile", help="BED file containing the amplicon scheme")
    parser.add_argument(
        "--normalise", type=int, help="Subsample to n coverage per strand"
    )
    parser.add_argument(
        "--min-mapq", type=int, default=20, help="Minimum mapping quality to keep"
    )
    parser.add_argument(
        "--primer-match-threshold",
        type=int,
        default=5,
        help="Fuzzy match primer positions within this threshold",
    )
    parser.add_argument("--report", type=str, help="Output report to file")
    parser.add_argument(
        "--amp-depth-report", type=str, help="Output amplicon depth tsv to path"
    )
    parser.add_argument(
        "--trim-primers",
        action="store_true",
        help="Trims primers from reads",
    )
    parser.add_argument("--paired", action="store_true", help="Paired end reads")
    parser.add_argument(
        "--no-read-groups",
        dest="no_read_groups",
        action="store_true",
        help="Do not divide reads into groups in SAM output",
    )
    parser.add_argument("--verbose", action="store_true", help="Debug mode")
    parser.add_argument("--remove-incorrect-pairs", action="store_true")

    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()
