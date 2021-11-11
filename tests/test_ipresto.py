import os
import subprocess
# import pytest


def test_ipresto_help(tmp_path):
    """Tests --help from main ipresto.py script from command line"""
    e = subprocess.check_call("python ipresto.py -h")
    assert e == 0, "Help message for ipresto.py failed"


def test_ipresto(tmp_path):
    """Tests main ipresto.py script from command line"""
    pass


if __name__ == '__main__':
    test_dir = './'
    test_ipresto_help(test_dir)
    test_ipresto(test_dir)
