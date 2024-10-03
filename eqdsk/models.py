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
        """
        Returns
        -------
        :
            The difference of the values as a new ZeroOne object

        Raises
        ------
        TypeError
            if not a ZeroOne object
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

    @classmethod
    def _missing_(cls, value) -> Sign:
        """

        Returns
        -------
        :
            The Sign object corresponding to the value.

        Raises
        ------
        ValueError
            If the value is not a known Sign.
        """
        if isinstance(value, bool):
            return Sign.POSITIVE if value else Sign.NEGATIVE
        raise ValueError(f"'{value}' not a known Sign") from None

    def __mul__(self, other: Any) -> Sign | int:
        """
        Returns
        -------
        :
            The product of the sign with the other value.

        Notes
        -----
        - If it is another Sign, return the product of the values.
        - If it is a number, return the product of the value and the number.
        """
        if type(other) is Sign:
            return Sign(self.value * other.value)
        return self.value * other
