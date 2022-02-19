import filecmp
import os
import shutil
import tempfile
import pathlib
import difflib
import webbrowser
import subprocess

import pytest
from collections import namedtuple


@pytest.mark.integration_tests
def test_eg01_rplot():

    test_a01 = """Rscript source Scripts/SwitchErrorTest_PhaseExtenderSetA_02.R"""

    subprocess.run(test_a01, shell=True)

    assert (
        is_same_file(
            "tests/testfiles/TestFullApp_output/viewvcfout01.dict",
            "tests/testfiles/TestFullApp_reference/ref_viewvcfout01.dict",
        )
        is True
    )

