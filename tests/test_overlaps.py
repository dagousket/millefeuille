import os
import sys
import pytest
import shutil
import pyranges as pr

from millefeuille.module import overlaps as ov


def test_load_beds():
    # Test with three bed files
    file_beds = ["./tests/sample1.bed", "./tests/sample2.bed", "./tests/sample3.bed"]
    result = ov.load_beds(file_beds)
    assert len(result) == 3, "Expected 3 bed files to be loaded"
    assert all(
        isinstance(bed, pr.PyRanges) for bed in result
    ), "Each list element should be a PyRanges object"
