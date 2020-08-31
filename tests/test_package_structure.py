"""ensure that canela package and subpackages can be correctly imported"""
import pytest


def test_package():
    import canela
    assert isinstance(canela.__version__, str)


def test_bin():
    import canela.bin

