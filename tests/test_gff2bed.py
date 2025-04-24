import os
import sys
import pytest
import shutil

from millefeuille.module import gff2bed as g2b

def test_get_featureDict():
  gff_entry = "gene_id=FBgn0002121;parent_type=mRNA;Name=l(2)gl-cds;Parent=FBtr0306591,FBtr0330655"
  expected_output = {'gene_id': ['FBgn0002121'], 'parent_type': ['mRNA'], 'Name': ['l(2)gl-cds'], 'Parent': ['FBtr0306591', 'FBtr0330655']}
  result = g2b.get_featureDict(gff_entry)
  assert isinstance(result, dict)
  assert result == expected_output, f"Expected {expected_output}, but got {result}"

def test_get_Dictgff():
  file_gff = "./tests/sample.gff"
  result = g2b.get_Dictgff(file_gff)
  expected_output = {
    'FBtr0078049': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}],
    'FBtr0078050': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}],
    'FBtr0078047': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}],
    'FBtr0078048': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}],
    'FBtr0330674': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}]
    }
  assert isinstance(result, dict)
  assert result == expected_output, f"Expected {expected_output}, but got {result}"

def test_consistency_check_ok(capsys):
  god_gff = {'FBtr0078049': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}], 'FBtr0078050': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}], 'FBtr0078047': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}], 'FBtr0078048': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}], 'FBtr0330674': [{'chr': 'chr2L', 'start': '337016', 'stop': '337198', 'strand': '+', 'name': 'gene_id=FBgn0004611;parent_type=mRNA;Name=Plc21C:9;Parent=FBtr0078049,FBtr0078050,FBtr0078047,FBtr0078048,FBtr0330674'}]}
  result = g2b.consistency_check(god_gff)
  expected_output = True
  captured = capsys.readouterr()
  assert "Discarded 0 elements" in captured.out
  assert isinstance(result, bool)
  assert result == expected_output, f"Expected {expected_output}, but got {result}"


def test_consistency_check_nok(capsys):
  bad_gff = {
    'FBtr0084082': [{'chr': '3R', 'start': '21367910', 'stop': '21368238', 'strand': '+', 'name': 'Parent=FBtr0084082'}, {'chr': '3R', 'start': '21375060', 'stop': '21375912', 'strand': '-', 'name': 'Parent=FBtr0084066'}],
    'FBtr0084066': [{'chr': '3R', 'start': '21376602', 'stop': '21376900', 'strand': '-', 'name': 'Parent=FBtr0084066'}, {'chr': '3R', 'start': '21376819', 'stop': '21377076', 'strand': '-', 'name': 'Parent=FBtr0084066'}],
    'FBtr0084067': [{'chr': '3L', 'start': '21375060', 'stop': '21375912', 'strand': '-', 'name': 'Parent=FBtr0084066'}, {'chr': '3R', 'start': '21376602', 'stop': '21376741', 'strand': '-', 'name': 'Parent=FBtr0084066'}]}
  result = g2b.consistency_check(bad_gff)
  captured = capsys.readouterr()
  assert "Discarded 2 elements" in captured.out
  assert "1 strand inconsistencies" in captured.out
  assert "1 Parent with overlapping exon" in captured.out
  assert "1 chromosome inconsistencies" in captured.out
  modified_gff = {'FBtr0084067': [{'chr': '3L', 'start': '21375060', 'stop': '21375912', 'strand': '-', 'name': 'Parent=FBtr0084066'}, {'chr': '3R', 'start': '21376602', 'stop': '21376741', 'strand': '-', 'name': 'Parent=FBtr0084066'}]}
  assert isinstance(bad_gff, dict)
  assert bad_gff == modified_gff, f"Expected {modified_gff}, but got {bad_gff}"

def test_bed6_generator(tmpdir) :
  file_gff = os.path.abspath("./tests/sample.gff")
  with tmpdir.as_cwd() as old_dir:
      shutil.copy(file_gff, os.path.join(os.getcwd(), "sample.gff"))
      result = g2b.bed6_generator(file_gff = "sample.gff", bedname = "sample")
      # Check if the GFF file was created
      bed_file = os.path.join(os.getcwd(), "sample.bed6")
      assert os.path.exists(bed_file), f"BED file {bed_file} was not created."
      with open(bed_file, 'r') as f:
          lines = f.readlines()
          assert len(lines) > 0, "BED file is empty."
      # Check if the BED file has the correct number of columns
      for line in lines:
          columns = line.strip().split('\t')
          assert len(columns) == 6, f"BED file {bed_file} does not have 6 columns."
      # Check if the BED file has the correct format
      for line in lines:
          columns = line.strip().split('\t')
          assert columns[0].startswith("chr"), f"BED file {bed_file} does not have the correct format."
          assert columns[1].isdigit(), f"BED file {bed_file} does not have the correct format."
          assert columns[2].isdigit(), f"BED file {bed_file} does not have the correct format."
          assert columns[3].isascii(), f"BED file {bed_file} does not have the correct format."
          assert columns[5] in ['+', '-'], f"BED file {bed_file} does not have the correct format."
      # Check if the BED file has the correct number of lines
      assert len(lines) == 5, f"BED file {bed_file} does not have the correct number of lines."

def test_bed12_generator(tmpdir) :
  file_gff = os.path.abspath("./tests/sample.gff")
  with tmpdir.as_cwd() as old_dir:
      shutil.copy(file_gff, os.path.join(os.getcwd(), "sample.gff"))
      result = g2b.bed12_generator(file_gff = "sample.gff", bedname = "sample")
      # Check if the GFF file was created
      bed_file = os.path.join(os.getcwd(), "sample.bed12")
      assert os.path.exists(bed_file), f"BED file {bed_file} was not created."
      with open(bed_file, 'r') as f:
          lines = f.readlines()
          assert len(lines) > 0, "BED file is empty."
      # Check if the BED file has the correct number of columns
      for line in lines:
          columns = line.strip().split('\t')
          assert len(columns) == 12, f"BED file {bed_file} does not have 12 columns."
      # Check if the BED file has the correct format
      for line in lines:
          columns = line.strip().split('\t')
          assert columns[0].startswith("chr"), f"BED file {bed_file} does not have the correct format."
          assert columns[1].isdigit(), f"BED file {bed_file} does not have the correct format."
          assert columns[2].isdigit(), f"BED file {bed_file} does not have the correct format."
          assert columns[3].isascii(), f"BED file {bed_file} does not have the correct format."
          assert columns[5] in ['+', '-'], f"BED file {bed_file} does not have the correct format."
