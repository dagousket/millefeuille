import os
import sys
import argparse
import pyranges as pr, pandas as pd


def load_beds(list_bed: list) -> list:
    """
    Create a list of pyranges intervals from a list of bed files.

    Parameters
    ----------
    list_bed : list
      A list of bed files. Must be of length 3.

    Returns
    -------
    list
      A list of pyranges intervals.
    """
    if len(list_bed) != 3:
        raise ValueError("list_bed must be of length 3")
    if not all(os.path.isfile(bed) for bed in list_bed):
        raise ValueError("All elements of list_bed must be valid file paths")

    pr_list = [pr.readers.read_bed(bed) for bed in list_bed]

    return pr_list
