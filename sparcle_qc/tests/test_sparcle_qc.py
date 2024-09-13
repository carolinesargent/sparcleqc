"""
Unit and regression test for the sparcle_qc package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import sparcle_qc


def test_sparcle_qc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "sparcle_qc" in sys.modules
