import os
import subprocess
# import pytest


def test_ipresto_help(tmp_path):
    """Tests --help from main ipresto.py script from command line"""
    tests_dir = os.path.split(os.path.realpath(__file__))[0]
    ipresto_path = os.path.join(tests_dir, "../ipresto.py")
    e = subprocess.check_call(f"python {ipresto_path} -h", shell=True)
    assert e == 0, "Help message for ipresto.py failed"


def test_ipresto(tmp_path):
    """Tests main ipresto.py script from command line"""
    tests_dir = os.path.split(os.path.realpath(__file__))[0]
    ipresto_path = os.path.join(tests_dir, "../ipresto.py")
    test_files = os.path.join(tests_dir, "test_files")
    test_out = os.path.join(tmp_path, "out")
    hmm_path = os.path.join(
        tests_dir, "../../../subcluster_data/domains/Pfam_100subs_tc.hmm")
    include_list = os.path.join(tests_dir, "../files/biosynthetic_domains.txt")
    cmd = f"python {ipresto_path} -i {test_files} -o {test_out} " \
          f"--hmm_path {hmm_path} -c 4 --no_redundancy_filtering " \
          f"--include_list {include_list} --remove_genes_below_count 0 -p 1.0"
    print(cmd)
    e = subprocess.check_call(cmd, shell=True)
    assert e == 0, "Basic iPRESTO run failed on test data"
    # todo: add test for using include_list, remove_genes_below_count, etc


def test_ipresto_red_filtering(tmp_path):
    """Tests main ipresto.py script from command line"""
    tests_dir = os.path.dirname(os.path.realpath(__file__))
    ipresto_path = os.path.join(tests_dir, "../ipresto.py")
    test_files = os.path.join(tests_dir, "test_files")
    test_out = os.path.join(tmp_path, "out_red")
    hmm_path = os.path.join(
        tests_dir, "../../../subcluster_data/domains/Pfam_100subs_tc.hmm")
    include_list = os.path.join(tests_dir, "../files/biosynthetic_domains.txt")
    cmd = f"python {ipresto_path} -i {test_files} -o {test_out} " \
          f"--hmm_path {hmm_path} -c 4 " \
          f"--include_list {include_list} --remove_genes_below_count 0 -p 1.0"
    print(cmd)
    e = subprocess.check_call(cmd, shell=True)
    assert e == 0, "Basic iPRESTO run failed on test data"


if __name__ == '__main__':
    # use the test directory to make output files in (in out)
    test_out_dir = os.path.dirname(os.path.realpath(__file__))
    test_ipresto_help(test_out_dir)
    test_ipresto(test_out_dir)
    test_ipresto_red_filtering(test_out_dir)
