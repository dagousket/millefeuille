import os
import sys
import pytest

from millefeuille.module import bed2gff as b2g

def test_get_Dictbed():
    # Test with a sample BED file
    expected_output = [
      {'chr': 'chr1', 'start': '20', 'end': '40', 'features': 'peak_1', 'strand': '+'},
      {'chr': 'chr1', 'start': '60', 'end': '110', 'features': 'peak_2', 'strand': '+'},
      {'chr': 'chr1', 'start': '170', 'end': '180', 'features': 'peak_3', 'strand': '+'}
      ]
    file_bed = "./tests/sample1.bed"
    result = b2g.get_Dictbed(file_bed)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

def test_get_Dictbed12():
    # Test with a sample BED12 file
    expected_output = [
      {'chr': 'chr1', 'start': '20', 'features': 'peak_1', 'strand': '+', 'block': '1', 'size': '182', 'starting_block': '0'},
      {'chr': 'chr1', 'start': '60', 'features': 'peak_2', 'strand': '+', 'block': '1', 'size': '182', 'starting_block': '0'},
      {'chr': 'chr1', 'start': '170', 'features': 'peak_3', 'strand': '+', 'block': '1', 'size': '182', 'starting_block': '0'}
      ]
    file_bed12 = "./tests/sample1.bed12"
    result = b2g.get_Dictbed12(file_bed12)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"
