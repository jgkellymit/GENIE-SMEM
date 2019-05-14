from ExactMatch import ExactMatch
from LUT import LUT
from SMEM import SMEM
import datetime


# LUT tests
def test_full():
    match = ExactMatch("mississippi.fa")
    match.load_fm_index()

    q = "missi"

    smem = SMEM(match)
    print(smem.get_smems_lut(q))


def test_simple():
    match = ExactMatch("mississippi.fa")
    match.load_fm_index()

    q = "mmissi"

    smem = SMEM(match)
    print(smem.get_smems_lut(q))


# TODO breaks it
def test_none():
    match = ExactMatch("mississippi.fa")
    match.load_fm_index()

    q = "aaaaa"

    smem = SMEM(match)
    print(smem.get_smems_lut(q))


def test_case1():
    match = ExactMatch("mississippi.fa")
    match.load_fm_index()

    q = "missippi"

    smem = SMEM(match)
    print(smem.get_smems_lut(q))


# TODO
def test_case2():
    match = ExactMatch("mississippi.fa")
    match.load_fm_index()

    q = "missippi"

    smem = SMEM(match)
    print(smem.get_smems_lut(q))


def test_case3():
    match = ExactMatch("mississippi.fa")
    match.load_fm_index()

    q = "mmiss"

    smem = SMEM(match)
    print(smem.get_smems_lut(q))


def test_case4():
    match = ExactMatch("mississippi.fa")
    match.load_fm_index()

    q = "mmissippss"

    smem = SMEM(match)
    print(smem.get_smems_lut(q))


class foo(object):
    def __init__(self):
        pass

    def __len__(self):
        return -1

if __name__ == '__main__':
    # test_case4()
    a = ["qqq", "asda", "asdaeaf"]
    print(max(a))


