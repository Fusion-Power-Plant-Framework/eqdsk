# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Eqdsk logging"""

import logging

logger = logging.getLogger("EQDSK Logger")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)8s] %(message)s",
)


def eqdsk_warn(*args, **kwargs):
    """Warning"""
    stacklevel = kwargs.pop("stacklevel", 1) + 1
    logger.warning(*args, stacklevel=stacklevel, **kwargs)


def eqdsk_print(*args, **kwargs):
    """Printing"""
    stacklevel = kwargs.pop("stacklevel", 1) + 1
    logger.info(*args, stacklevel=stacklevel, **kwargs)


def eqdsk_banner():
    """
    Returns
    -------
    :
        The eqdsk banner
    """
    from eqdsk._version import version  # noqa: PLC0415

    return f"""
********************************
*    ___ ___  ___  ___ _  __   *
*   | __/ _ \\|   \\/ __| |/ /   *
*   | _| (_) | |) \\__ \\ ' <    *
*   |___\\__\\_\\___/|___/_|\\_\\   *
*                              *
*   The Bluemira Developers    *
********************************
EQDSK Version: {version}
"""
