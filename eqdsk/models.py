# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""Enum models for values used in the COCOS specification."""

from __future__ import annotations

from enum import Enum
from typing import Any


class ZeroOne(Enum):
    """An enum representing the values 0 and 1 for the 2pi exponent of Bp."""

    ZERO = 0
    ONE = 1

    def __sub__(self, other: Any) -> ZeroOne:
        """Return the difference between the value and the other value.

        - If it is another ZeroOne, return the difference of the values.
        - Raise a `TypeError` otherwise.
        """
        if type(other) is ZeroOne:
            return ZeroOne(self.value - other.value)
        raise TypeError(
            f"Cannot subtract {type(other)} from {type(self)}.",
        )


class Sign(Enum):
    """An enum representing the positive or negative sign of
    a COCOS parameter.
    """

    POSITIVE = 1
    NEGATIVE = -1

    def __mul__(self, other: Any):
        """Return the product of the sign with the other value.

        - If it is another Sign, return the product of the values.
        - If it is a number, return the product of the value and the number.
        """
        if type(other) is Sign:
            return Sign(self.value * other.value)
        return self.value * other
