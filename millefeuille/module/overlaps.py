import os
import sys
import argparse
import pyranges as pr, pandas as pd
from itertools import combinations


def load_beds(list_bed: list, names: list = ["a", "b", "c"]) -> dict:
    """
    Create a list of pyranges intervals from a list of bed files.

    Parameters
    ----------
    list_bed : list
      A list of bed files. Must be of length 3.

    Returns
    -------
    list
      A dict of pyranges intervals.
    """
    if len(list_bed) != 3:
        raise ValueError("list_bed must be of length 3")
    if not all(os.path.isfile(bed) for bed in list_bed):
        raise ValueError("All elements of list_bed must be valid file paths")

    pr_dict = dict(zip(names, [pr.readers.read_bed(bed) for bed in list_bed]))

    return pr_dict


def double_overlap(pr_dict: dict) -> dict:
    """
    Calculate the overlaps between the two pyranges intervals.
    Only returns regions that do not overlap with the third interval.

    Parameters
    ----------
    pr_dict : dict
      A dict of pyranges intervals.
    output : dict
      The dictionnary of intervals that overlap with the two pyranges intervals.
    """
    bedpair_output = dict()
    for bedpair in combinations(pr_dict.keys(), 2):
        pair_overlap = pr_dict[bedpair[0]].intersect(pr_dict[bedpair[1]])
        last_overlap = [bed for bed in pr_dict.keys() if bed not in bedpair][0]
        pair_overlap = pair_overlap.overlap(pr_dict[last_overlap], invert=True)
        bedpair_output["::".join(bedpair)] = pair_overlap

    return bedpair_output


def triple_overlap(pr_dict: dict) -> dict:
    """
    Calculate the overlaps between the three pyranges intervals.

    Parameters
    ----------
    pr_dict : dict
      A dict of pyranges intervals.
    output : dict
      The intervals that overlap with all three intervals.
    """
    pr_dict_val = list(pr_dict.values())
    triple_overlap = pr_dict_val[0].intersect(pr_dict_val[1]).intersect(pr_dict_val[2])
    output = {"::".join(pr_dict.keys()): triple_overlap}

    return output


def single_overlap(pr_dict: dict) -> dict:
    """
    Calculate the intervals that do not overlap any of the two other pyranges intervals.

    Parameters
    ----------
    pr_dict : dict
      A dict of pyranges intervals.
    output : dict
      The intervals that overlap with any of the two other intervals.
    """
    bedsingle_output = dict()
    for bedpair in combinations(pr_dict.keys(), 2):
        last_overlap = [bed for bed in pr_dict.keys() if bed not in bedpair][0]
        single_overlap = pr_dict[last_overlap].overlap(pr_dict[bedpair[0]], invert=True)
        single_overlap = single_overlap.overlap(pr_dict[bedpair[1]], invert=True)
        bedsingle_output[last_overlap] = single_overlap

    return bedsingle_output
