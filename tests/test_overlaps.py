import os
import sys
import pytest
import shutil
import pyranges as pr
import pandas as pd

from millefeuille.module import overlaps as ov


def test_load_beds():
    # Test with three bed files
    file_beds = ["./tests/sample1.bed", "./tests/sample2.bed", "./tests/sample3.bed"]
    result = ov.load_beds(file_beds)
    assert len(result) == 3, "Expected 3 bed files to be loaded"
    assert all(
        isinstance(bed, pr.PyRanges) for bed in result.values()
    ), "Each list element should be a PyRanges object"


def test_triple_overlap():
    file_beds = ["./tests/sample1.bed", "./tests/sample2.bed", "./tests/sample3.bed"]
    pr_dict = ov.load_beds(file_beds)
    result = ov.triple_overlap(pr_dict)
    assert len(result) == 1, "Expected 1 overlap region"
    assert isinstance(result, dict), "Expected result to be a dictionary"
    assert all(
        isinstance(bed, pr.PyRanges) for bed in result.values()
    ), "Each list element should be a PyRanges object"
    assert list(result.keys()) == ["a::b::c"], "Expected key to be 'a::b::c'"
    expect_overlap = pd.DataFrame(
        {
            "Chromosome": {0: "chr1", 1: "chr1"},
            "Start": {0: 30, 1: 70},
            "End": {0: 35, 1: 80},
            "Name": {0: "peak_1", 1: "peak_2"},
            "Score": {0: 1000, 1: 1000},
            "Strand": {0: "+", 1: "+"},
        }
    )
    assert all(
        result["a::b::c"].as_df() == expect_overlap
    ), "Expected overlap regions coordinates to match"


def test_double_overlap():
    file_beds = ["./tests/sample1.bed", "./tests/sample2.bed", "./tests/sample3.bed"]
    pr_dict = ov.load_beds(file_beds)
    result = ov.double_overlap(pr_dict)
    assert len(result) == 3, "Expected 3 overlap region pairs"
    assert isinstance(result, dict), "Expected result to be a dictionary"
    assert all(
        isinstance(bed, pr.PyRanges) for bed in result.values()
    ), "Each list element should be a PyRanges object"
    assert list(result.keys()) == [
        "a::b",
        "a::c",
        "b::c",
    ], "Expected key to be 'a::b' 'a::c', 'b::c'"
    expect_overlap_ab = pd.DataFrame(
        {
            "Chromosome": {0: "chr1"},
            "Start": {0: 100},
            "End": {0: 110},
            "Name": {0: "peak_2"},
            "Score": {0: 1000},
            "Strand": {0: "+"},
        }
    )
    assert all(
        result["a::b"].as_df() == expect_overlap_ab
    ), "Expected overlap regions coordinates to match between a and b"
    expect_overlap_bc = pd.DataFrame(
        {
            "Chromosome": {0: "chr1"},
            "Start": {0: 150},
            "End": {0: 160},
            "Name": {0: "peak_4"},
        }
    )
    assert all(
        result["b::c"].as_df() == expect_overlap_bc
    ), "Expected overlap regions coordinates to match between b and c"
    assert (
        len(result["a::c"]) == 0
    ), "Expected no overlap regions coordinates between a and c"


def test_single_overlap():
    file_beds = ["./tests/sample1.bed", "./tests/sample2.bed", "./tests/sample3.bed"]
    pr_dict = ov.load_beds(file_beds)
    result = ov.single_overlap(pr_dict)
    assert len(result) == 3, "Expected 3 overlap regions"
    assert isinstance(result, dict), "Expected result to be a dictionary"
    assert all(
        isinstance(bed, pr.PyRanges) for bed in result.values()
    ), "Each list element should be a PyRanges object"
    assert set(result.keys()) == set(
        [
            "a",
            "b",
            "c",
        ]
    ), "Expected key to be 'a', 'b', 'c'"
    expect_overlap_a = pd.DataFrame(
        {
            "Chromosome": {0: "chr1"},
            "Start": {0: 170},
            "End": {0: 180},
            "Name": {0: "peak_3"},
            "Score": {0: 1000},
            "Strand": {0: "+"},
        }
    )
    assert all(
        result["a"].as_df() == expect_overlap_a
    ), "Expected overlap regions coordinates to match for a"
    assert len(result["b"]) == 0, "Expected no overlap regions coordinates for b"
    assert len(result["c"]) == 0, "Expected no overlap regions coordinates for c"


def test_all_overlaps():
    file_beds = ["./tests/sample1.bed", "./tests/sample2.bed", "./tests/sample3.bed"]
    result = ov.all_overlaps(file_beds)
    assert len(result) == 7, "Expected 7 overlap regions"
    assert isinstance(result, dict), "Expected result to be a dictionary"
    assert set(result.keys()) == set(
        [
            "a::b",
            "a::c",
            "b::c",
            "a::b::c",
            "a",
            "b",
            "c",
        ]
    ), "Expected key to be 'a::b', 'a::c', 'b::c', 'a::b::c', 'a', 'b', 'c'"
    assert result == {
        "c": 0,
        "b": 0,
        "a": 1,
        "a::b": 1,
        "a::c": 0,
        "b::c": 1,
        "a::b::c": 2,
    }, "Expect correct detection of overlaps"
