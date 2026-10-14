"""The simplest possible test: can Python find and load the package?"""

import signatomics
import signatomics.olink


def test_package_has_a_version():
    assert signatomics.__version__ == "0.0.1"
