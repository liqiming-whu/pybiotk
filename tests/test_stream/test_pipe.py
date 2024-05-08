import pytest
from stream.pipe import *


@pytest.fixture(scope="module")
def test_case1():
    return [1, 2, 3, 4, 5]

@pytest.fixture(scope="module")
def test_case2():
    return [1, 2, 2, 3, 4, 5, 1]

def test_head(test_case1):
    assert test_case1 | head(n=1) | to_list == [1]
    

def test_tail(test_case1):
    assert test_case1 | tail(n=2) | to_list == [4, 5]
    
    
def test_skip(test_case1):
    assert test_case1 | skip(n=2) | to_list == [3, 4, 5]
    

def test_filter(test_case1):
    assert test_case1 | filter(lambda x: x % 2 == 0) | to_list == [2, 4]


def test_dedup(test_case2):
    assert test_case2 | dedup(key=lambda x: x) | to_list == [1, 2, 3, 4, 5]

def test_uniq(test_case2):
    assert test_case2 | uniq | to_list == [1, 2, 3, 4, 5, 1]

def test_count(test_case2):
    assert test_case2 | count == 7

def test_permutations():
    assert "ABCD" | permutations(r=2) | mapwith("".join) | to_list == ["AB", "AC", "AD", "BA", "BC", "BD", "CA", "CB", "CD", "DA", "DB", "DC"]
    assert range(3) | permutations | mapwith(lambda x: "".join(map(str, x))) | to_list == ["012", "021", "102", "120", "201", "210"]
    
def test_join():
    assert [1, 2, 3] | join("|") == "1|2|3"
    assert ["1", "2", "3"] | join("|") == "1|2|3"

def test_flatten():
    assert [[1, 2], [3, 4]] | flatten | to_list== [1, 2, 3, 4]
    assert [1, 2, 3, 4] | flatten | to_list == [1, 2, 3, 4]
    assert [[1, 2], [[3, 4], [5, 6]]] | flatten | to_list == [1, 2, 3, 4, 5, 6]