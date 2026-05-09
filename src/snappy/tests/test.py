from snappy import testing
import snappy.tests
import snappy.tests.orb

modules = [
    snappy.tests.orb
]

def run_doctests(verbose=False, print_info=True):
    return testing.doctest_modules(modules,
                                   verbose=verbose,
                                   print_info=print_info)
run_doctests.__name__ = snappy.tests.__name__

if __name__ == '__main__':
    testing.run_doctests_as_main(run_doctests)
