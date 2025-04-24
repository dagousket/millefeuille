import os
import sys
import pytest
import shutil

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

def test_get_gff(tmpdir):
  file_bed = os.path.abspath("./tests/sample1.bed")
  with tmpdir.as_cwd() as old_dir:
      shutil.copy(file_bed, os.path.join(os.getcwd(), "sample1.bed"))
      result = b2g.get_gff(file_bed = "sample1.bed",
                       source = "test_source",
                       mol_type = "test_mol_type",
                       make_gff3 = True)
      # Check if the GFF file was created
      gff_file = os.path.join(os.getcwd(), "sample1.gff")
      assert os.path.exists(gff_file), f"GFF file {gff_file} was not created."
      # Check the content of the GFF file
      with open(gff_file, 'r') as f:
          lines = f.readlines()
          assert len(lines) > 0, "GFF file is empty."
          # Check the first line
          first_line = lines[0].strip().split('\t')
          assert len(first_line) == 9, "First line of GFF file does not have 9 columns."
          assert first_line[0] == "chr1", "First column of GFF file is not 'chr1'."
          assert first_line[2] == "test_mol_type", "Third column of GFF file is not 'test_mol_type'."
          assert first_line[6] == "+", "Seventh column of GFF file is not '+'."
          assert first_line[8] == "feature_id=peak_1", "Ninth column of GFF file is not 'feature_id=peak_1'."

def test_get_gff_from_bed12(tmpdir):
  file_bed = os.path.abspath("./tests/sample1.bed12")
  with tmpdir.as_cwd() as old_dir:
      shutil.copy(file_bed, os.path.join(os.getcwd(), "sample1.bed12"))
      result = b2g.get_gff(file_bed = "sample1.bed12",
                       source = "test_source",
                       mol_type = "test_mol_type",
                       make_gff3 = True)
      # Check if the GFF file was created
      gff_file = os.path.join(os.getcwd(), "sample1.gff")
      assert os.path.exists(gff_file), f"GFF file {gff_file} was not created."
      with open(gff_file, 'r') as f:
          lines = f.readlines()
          assert len(lines) > 0, "GFF file is empty."
          # Check the first line
          first_line = lines[0].strip().split('\t')
          assert len(first_line) == 9, "First line of GFF file does not have 9 columns."
          assert first_line[0] == "chr1", "First column of GFF file is not 'chr1'."
          assert first_line[2] == "test_mol_type", "Third column of GFF file is not 'test_mol_type'."
          assert first_line[6] == "+", "Seventh column of GFF file is not '+'."
          assert first_line[8] == "feature_id=peak_1", "Ninth column of GFF file is not 'feature_id=peak_1'."
