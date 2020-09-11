"""ensure that canela package and subpackages can be correctly imported"""
import pytest


def test_package():
    import canela
    assert isinstance(canela.__version__, str)


def test_subpackages():
    import canela.bin
    import canela.io
    import canela.data
    import canela.lpnc

