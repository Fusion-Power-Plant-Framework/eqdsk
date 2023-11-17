# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""Definitions for the COCOS system."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from enum import Enum, unique
from typing import TYPE_CHECKING, Any, Optional

import numpy as np

if TYPE_CHECKING:
    from eqdsk.file import EQDSKInterface


class ZeroOne(Enum):
    """An enum representing the values 0 and 1 of the exponent for Bp."""

    ZERO = 0
    ONE = 1

    def __sub__(self, other: Any) -> ZeroOne:  # noqa: ANN401
        if type(other) is ZeroOne:
            return ZeroOne(self.value - other.value)
        raise TypeError(
            f"Cannot subtract {type(other)} from {type(self)}.",
        )


class Sign(Enum):
    """
    An enum representing the
    positive or negative sign for a COCOS parameter.
    """

    POSITIVE = 1
    NEGATIVE = -1

    def __mul__(self, other: Any):  # noqa: ANN401
        if type(other) is Sign:
            return Sign(self.value * other.value)
        return self.value * other


@dataclass(frozen=True)
class COCOSParams:
    """The parameters for a single COCOS definition."""

    index: int
    exp_Bp: ZeroOne
    sign_Bp: Sign
    sign_R_phi_Z: Sign
    sign_rho_theta_phi: Sign


@unique
class COCOS(Enum):
    """An enum representing the 16 COCOS definitions."""

    C1 = COCOSParams(
        index=1,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C11 = COCOSParams(
        index=11,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C2 = COCOSParams(
        index=2,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C12 = COCOSParams(
        index=12,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C3 = COCOSParams(
        index=3,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C13 = COCOSParams(
        index=13,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C4 = COCOSParams(
        index=4,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C14 = COCOSParams(
        index=14,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C5 = COCOSParams(
        index=5,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C15 = COCOSParams(
        index=15,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C6 = COCOSParams(
        index=6,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C16 = COCOSParams(
        index=16,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C7 = COCOSParams(
        index=7,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C17 = COCOSParams(
        index=17,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C8 = COCOSParams(
        index=8,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C18 = COCOSParams(
        index=18,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )

    def __init__(self, c: COCOSParams):
        # shortcuts .value access
        self.index = c.index
        self.exp_Bp = c.exp_Bp
        self.sign_Bp = c.sign_Bp
        self.sign_R_phi_Z = c.sign_R_phi_Z
        self.sign_rho_theta_phi = c.sign_rho_theta_phi

    @classmethod
    def with_index(cls, cocos_index: int) -> COCOS:
        """Return the COCOS with the given index."""
        if not (cocos_index in range(1, 9) or cocos_index in range(11, 19)):
            msg = f"Convention number {cocos_index} is not valid. "
            "Must be between 1 and 8 or 11 and 18."
            raise ValueError(msg)
        return next(x for x in cls if x.index == cocos_index)

    @classmethod
    def matching_convention(
        cls,
        exp_Bp: ZeroOne,
        sign_Bp: Sign,
        sign_R_phi_Z: Sign,
        sign_rho_theta_phi: Sign,
    ) -> COCOS:
        """Return the COCOS matching the given parameters."""

        def _match_cocos(c: COCOS) -> bool:
            return all(
                [
                    c.exp_Bp == exp_Bp,
                    c.sign_Bp == sign_Bp,
                    c.sign_R_phi_Z == sign_R_phi_Z,
                    c.sign_rho_theta_phi == sign_rho_theta_phi,
                ],
            )

        return next(filter(_match_cocos, cls))


@dataclass(frozen=True)
class COCOSTransform:
    """The values for a transformation between in and out COCOS definitions."""

    in_cocos: COCOS
    out_cocos: COCOS

    plasma_current: float
    b_toroidal: float
    coil_current: float
    poloidal_current: float
    psi: float
    pprime: float
    ffprime: float
    q: float


def identify_cocos(
    plasma_current: float,
    b_toroidal: float,
    psi_at_boundary: float,
    psi_at_mag_axis: float,
    q_psi: np.ndarray,
    *,
    phi_clockwise_from_top: bool,
    volt_seconds_per_radian: bool,
) -> COCOS:
    """
    Identify the COCOS for the given values from an eqdsk file.

    Parameters
    ----------
    plasma_current :
        The plasma current.
    b_toroidal :
        The toroidal magnetic field.
    psi_at_boundary :
        The psi value at the boundary.
    psi_at_mag_axis :
        The psi value at the magnetic axis.
    q_psi :
        The psi q (safety factor) values.
    phi_clockwise_from_top :
        Whether phi is clockwise from the top.
    volt_seconds_per_radian :
        Whether the flux is in volt seconds per radian.
    """
    sign_R_phi_Z = Sign.NEGATIVE if phi_clockwise_from_top else Sign.POSITIVE
    exp_Bp = ZeroOne.ONE if volt_seconds_per_radian else ZeroOne.ZERO

    sign_Ip = Sign(np.sign(plasma_current))
    sign_B0 = Sign(np.sign(b_toroidal))
    sign_psi_inc_towards_boundary = Sign(
        np.sign(psi_at_boundary - psi_at_mag_axis)
    )

    sign_q = np.sign(q_psi)
    if sign_q.min() != sign_q.max():
        raise ValueError(
            "The sign of qpsi is not consistent across the flux surfaces.",
        )
    sign_q = Sign(sign_q.max())

    sign_Bp = sign_Ip * sign_psi_inc_towards_boundary
    sign_rho_theta_phi = sign_Ip * sign_B0 * sign_q * sign_R_phi_Z

    return COCOS.matching_convention(
        exp_Bp=exp_Bp,
        sign_Bp=sign_Bp,
        sign_R_phi_Z=sign_R_phi_Z,
        sign_rho_theta_phi=sign_rho_theta_phi,
    )


def identify_eqdsk(
    eqdsk: EQDSKInterface,
    clockwise_phi: Optional[bool] = None,
    volt_seconds_per_radian: Optional[bool] = None,
) -> list[COCOS]:
    """
    Identify the COCOS for the given eqdsk file.

    Parameters
    ----------
    eqdsk :
        The eqdsk file to identify the COCOS for.
    clockwise_phi :
        Whether phi is clockwise from the top,
        by default None which means either.
    volt_seconds_per_radian :
        Whether the flux is in volt seconds per radian,
        by default None which means either.
    """
    if eqdsk.qpsi is None:
        raise ValueError("qpsi is not defined in the eqdsk file.")

    cw_phis = [True, False] if clockwise_phi is None else [clockwise_phi]
    vs_prs = (
        [True, False]
        if volt_seconds_per_radian is None
        else [volt_seconds_per_radian]
    )

    conventions = [
        identify_cocos(
            plasma_current=eqdsk.cplasma,
            b_toroidal=eqdsk.bcentre,
            psi_at_boundary=eqdsk.psibdry,
            psi_at_mag_axis=eqdsk.psimag,
            q_psi=eqdsk.qpsi,
            phi_clockwise_from_top=cw_phi,
            volt_seconds_per_radian=vs_pr,
        )
        for cw_phi in cw_phis
        for vs_pr in vs_prs
    ]

    # return sort asc by index
    conventions.sort(key=lambda x: x.index)
    return conventions


def transform_cocos(
    from_cocos_index: int, to_cocos_index: int
) -> COCOSTransform:
    """Return a transformation needed to transform from one COCOS to another."""
    in_cocos = COCOS.with_index(from_cocos_index)
    out_cocos = COCOS.with_index(to_cocos_index)

    eff_bp = (out_cocos.sign_Bp * in_cocos.sign_Bp).value
    eff_R_phi_Z = (out_cocos.sign_R_phi_Z * in_cocos.sign_R_phi_Z).value
    eff_rho_theta_phi = (
        out_cocos.sign_rho_theta_phi * in_cocos.sign_rho_theta_phi
    ).value
    eff_exp_bp = (out_cocos.exp_Bp - in_cocos.exp_Bp).value

    pi_factor = (2 * np.pi) ** eff_exp_bp

    return COCOSTransform(
        in_cocos=in_cocos,
        out_cocos=out_cocos,
        plasma_current=eff_R_phi_Z,
        b_toroidal=eff_R_phi_Z,
        coil_current=eff_R_phi_Z,
        poloidal_current=eff_rho_theta_phi,
        psi=eff_bp * eff_R_phi_Z * (1 / pi_factor),
        pprime=eff_bp * eff_R_phi_Z * pi_factor,
        ffprime=eff_bp * eff_R_phi_Z * pi_factor,
        q=eff_R_phi_Z * eff_rho_theta_phi,
    )


def convert_eqdsk(eqdsk: EQDSKInterface, to_cocos_index: int) -> EQDSKInterface:
    """Convert an eqdsk file to the given COCOS."""
    in_eqdsk = eqdsk
    out_eqdsk = deepcopy(in_eqdsk)

    in_index = eqdsk.cocos.index
    out_index = to_cocos_index

    if in_index == out_index:
        # returns the copy
        return out_eqdsk

    transform = transform_cocos(
        from_cocos_index=in_index,
        to_cocos_index=out_index,
    )

    out_eqdsk.cplasma = transform.plasma_current * in_eqdsk.cplasma
    out_eqdsk.bcentre = transform.b_toroidal * in_eqdsk.bcentre

    out_eqdsk.Ic = transform.coil_current * in_eqdsk.Ic
    out_eqdsk.fpol = transform.poloidal_current * in_eqdsk.fpol

    out_eqdsk.psi = transform.psi * in_eqdsk.psi
    out_eqdsk.psibdry = transform.psi * in_eqdsk.psibdry
    out_eqdsk.psimag = transform.psi * in_eqdsk.psimag

    out_eqdsk.pprime = transform.pprime * in_eqdsk.pprime
    out_eqdsk.ffprime = transform.ffprime * in_eqdsk.ffprime

    if in_eqdsk.qpsi is not None:
        out_eqdsk.qpsi = transform.q * in_eqdsk.qpsi

    # reidentify the eqdsk
    out_eqdsk.identify(
        clockwise_phi=transform.out_cocos.sign_R_phi_Z == Sign.NEGATIVE,
        volt_seconds_per_radian=transform.out_cocos.exp_Bp == ZeroOne.ONE,
    )
    if out_eqdsk.cocos.index != to_cocos_index:
        raise RuntimeError(
            f"Failed to convert eqdsk to COCOS {to_cocos_index}, "
            f"eqdsk file identifed as COCOS {out_eqdsk.cocos.index} "
            "post conversion.",
        )

    return out_eqdsk
