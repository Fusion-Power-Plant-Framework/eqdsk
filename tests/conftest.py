# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
from contextlib import suppress

import pytest

from tests._helpers import get_private_dir


def pytest_addoption(parser):
    """
    Adds a custom command line option to pytest.
    """
    parser.addoption(
        "--private",
        action="store_true",
        dest="private",
        default=False,
        help="run tests that use private data",
    )
    parser.addoption(
        "--plotting-on",
        action="store_true",
        default=False,
        help="switch on interactive plotting in tests",
    )


def pytest_configure(config):
    """
    Configures pytest.

    Raises
    ------
    ValueError
        Failed to find private data directory
    """
    if config.option.private and get_private_dir() is None:
        raise ValueError("You cannot run private tests. Data directory not found")

    config.option.markexpr = "" if config.option.private else "not private"

    if not config.option.plotting_on:
        # We're not displaying plots so use a display-less backend
        with suppress(ImportError):
            import matplotlib as mpl  # noqa: PLC0415

            mpl.use("Agg")


@pytest.fixture(autouse=True)
def _plot_show_and_close(request):
    """Fixture to show and close plots

    Notes
    -----
    Does not do anything if testclass marked with 'classplot'
    """
    with suppress(ImportError):
        import matplotlib.pyplot as plt  # noqa: PLC0415

        cls = request.node.getparent(pytest.Class)

        yield
        clstitle = "" if cls is None else cls.name
        for fig in list(map(plt.figure, plt.get_fignums())):
            fig.suptitle(
                f"{fig.get_suptitle()} {clstitle}::"
                f"{request.node.getparent(pytest.Function).name}"
            )
        plt.show()
        plt.close()
