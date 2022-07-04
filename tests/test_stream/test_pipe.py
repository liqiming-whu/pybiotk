import pytest
from stream.pipe import *


@pytest.fixture(scope="module")
def test_case():
    return [1, 2, 3, 4, 5]


def test_head(test_case):
    assert test_case | head(n=1) | to_list == [1]