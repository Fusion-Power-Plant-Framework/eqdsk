# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""The eqdsk package provides a class for reading and writing eqdsk files.

The core class is the [`EQDSKInterface`][eqdsk.file],
responsible for most of the functionality.

There is also the COCOS implementation in the [`cocos.py`][eqdsk.cocos] file,
which provides the coordinate system conversion functions.
"""

from eqdsk.file import EQDSKInterface
