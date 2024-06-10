# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

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


def pytest_configure(config):
    """
    Configures pytest.
    """
    if config.option.private and get_private_dir() is None:
        raise ValueError("You cannot run private tests. Data directory not found")

    config.option.markexpr = "" if config.option.private else "not private"
