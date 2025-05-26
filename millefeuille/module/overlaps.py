import os
import sys
import argparse
import pyranges as pr, pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
import upsetplot as upset
from matplotlib_venn import venn3


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


def all_overlaps(
    list_bed: list, names: list = ["a", "b", "c"], as_bp: bool = False
) -> dict:
    """
    Calculate the overlaps between the three pyranges intervals.

    Parameters
    ----------
    list_bed : list
      A list of bed files. Must be of length 3.
    names : list
      A list of names for the bed files. Must be of length 3.
    as_bp : bool
      If True, return the length of the intervals in base pairs instead of overlap count. Default is False.

    Returns
    -------
    list
      A dict of pyranges intervals.
    """
    all_beds = load_beds(list_bed)
    bed_1layer = single_overlap(all_beds)
    bed_2layer = double_overlap(all_beds)
    bed_3layer = triple_overlap(all_beds)
    all_overlap = {**bed_1layer, **bed_2layer, **bed_3layer}
    if as_bp:
        all_overlap = dict(
            zip(
                all_overlap.keys(),
                [ov.merge(strand=False).length for ov in all_overlap.values()],
            )
        )
    else:
        all_overlap = dict(
            zip(all_overlap.keys(), [len(ov) for ov in all_overlap.values()])
        )
    return all_overlap


def plot_overlaps(
    list_bed: list,
    names: list = ["a", "b", "c"],
    as_venn: bool = False,
    as_bp: bool = False,
) -> None:
    """
    Plot the overlaps between the three pyranges intervals.

    Parameters
    ----------
    list_bed : list
      A list of bed files. Must be of length 3.
    names : list
      A list of names for the bed files. Must be of length 3.
    as_venn : bool
      If True, plot a Venn diagram instead of an Upset plot. Default is False.
    as_bp : bool
      If True, return the length of the intervals in base pairs instead of overlap count. Default is False.
    ax : plt.Axes
      The matplotlib Axes to plot on. If None, a new figure and axes will be created. Default is None.
    """

    all_overlap = all_overlaps(list_bed, names, as_bp)

    if as_venn:
        # plot Venn with count order : Abc, aBc, ABc, abC, AbC, aBC, ABC
        ordered_items = [
            [x, y, [x, y], z, [x, z], [y, z], [x, y, z]]
            for x, y, z in zip(names[0], names[1], names[2])
        ][0]
        ordered_values = [all_overlap.get("::".join(item), 0) for item in ordered_items]
        venn3(subsets=ordered_values, set_labels=names)

    elif not as_venn:
        # plot Upset
        pd.set_option("future.no_silent_downcasting", True)
        counts = upset.from_memberships(
            [name.split("::") for name in all_overlap.keys()],
            data=list(all_overlap.values()),
        )
        upset.plot(counts)

    plt.show()

    return None
