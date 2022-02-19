import filecmp
import os
import shutil
import tempfile
import pathlib
import difflib
import webbrowser

import pytest
from collections import namedtuple




def test_first_example():
    with tempfile.TemporaryDirectory() as temp_dir1:
        soi_meta = soi_meta_cls('ms02g', [('MA605:PI', 'MA605:PG_al')], [('Sp21:PI', 'Sp21:PG_al')], 'pat_hap', 'mat_hap')
        args1 = ('tests/inputs/haplotype_file01.txt', soi_meta, temp_dir1, '', 1, 3, '+', 'no')
        phase_stich(*args1)
        assert is_same_dir('tests/outdir/egout1/', temp_dir1) 
