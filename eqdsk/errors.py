# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""Custom errors"""

from eqdsk.cocos import COCOS


class NoSingleConventionError(Exception):
    """Raised when no single COCOS can be determined from the EQDSK file."""

    def __init__(self, conventions: list[COCOS], message_extra: str):
        self.conventions = conventions
        self.message_extra = message_extra
        super().__init__(
            f"A single COCOS could not be determined, "
            f"found conventions "
            f"({', '.join([str(c.index) for c in self.conventions])}) "
            f"for the EQDSK file.\n{self.message_extra}"
        )


class MissingQpsiDataError(Exception):
    """Raised when attempting to identify the COCOS for an EQDSK file,
    that does not have an qpsi data in it.
    """

    def __init__(self, message_extra: str):
        self.message_extra = message_extra
        super().__init__(
            f"In order to properly identify the COCOS of this EQDSK file, "
            f"qpsi data must be present in the file.\n{self.message_extra}"
        )
