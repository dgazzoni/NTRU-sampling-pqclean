"""
Checks that no dynamic memory functions are used
"""

import pqclean
import helpers
import sys
import unittest


def test_dynamic_memory():
    if sys.platform not in ['linux', 'darwin']:
        raise unittest.SkipTest()
    for scheme in pqclean.Scheme.all_schemes():
        for implementation in scheme.implementations:
            # Keep this loop outside, to allow multiple assertions
            for function in ['malloc', 'free', 'realloc', 'calloc']:
                yield (check_dynamic_memory,
                       scheme.name, implementation.name, function)


def check_dynamic_memory(scheme_name, implementation_name, function):
    implementation = pqclean.Implementation.by_name(
        scheme_name, implementation_name)
    # 'make' will take care of not rebuilding existing library files
    helpers.run_subprocess(
        ['make'],
        implementation.path()
    )
    out = helpers.run_subprocess(
        ['nm', '-g', 'lib{}_{}.a'.format(scheme_name,
                                         implementation_name)],
        implementation.path()
    )

    lines = out.strip().split("\n")

    for line in lines:
        if 'U {}'.format(function) in line:
            raise AssertionError(
                "Illegal use of dynamic memory function '{}'".format(function))

if __name__ == '__main__':
    try:
        import nose2
        nose2.main()
    except ImportError:
        import nose
        nose.runmodule()