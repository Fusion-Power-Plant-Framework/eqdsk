from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from eqdsk.file import EQDSKInterface


class ZeroOne(Enum):
    ZERO = 0
    ONE = 1


class Sign(Enum):
    P = 1
    N = -1


@dataclass(frozen=True)
class COCOSConvention:
    cc_index: int
    exp_Bp: ZeroOne  # 0 -> ψ/2pi^(1-0) i.e. ψ/2pi, 1 -> ψ/2pi^(1-1) i.e. just ψ
    sign_Bp: Sign  # ψ (P -> increasing)/(N -> decreasing) from the magnetic axis
    sign_R_phi_Z: Sign  # P -> (R, ϕ, Z) ϕ cnt-clockwise, N -> (R, Z, ϕ) ϕ clockwise
    sign_rho_theta_phi: Sign  # P -> (ρ, θ, ϕ), N -> (ρ, ϕ, θ)


C1 = COCOSConvention(
    cc_index=1,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
)
C11 = COCOSConvention(
    cc_index=11,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
)
C2 = COCOSConvention(
    cc_index=2,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
)
C12 = COCOSConvention(
    cc_index=12,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
)
C3 = COCOSConvention(
    cc_index=3,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
)
C13 = COCOSConvention(
    cc_index=13,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
)
C4 = COCOSConvention(
    cc_index=4,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
)
C14 = COCOSConvention(
    cc_index=14,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
)
C5 = COCOSConvention(
    cc_index=5,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
)
C15 = COCOSConvention(
    cc_index=15,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
)
C6 = COCOSConvention(
    cc_index=6,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
)
C16 = COCOSConvention(
    cc_index=16,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
)
C7 = COCOSConvention(
    cc_index=7,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
)
C17 = COCOSConvention(
    cc_index=17,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
)
C8 = COCOSConvention(
    cc_index=8,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
)
C18 = COCOSConvention(
    cc_index=18,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
)


def get_convention(cocos_index: int) -> COCOSConvention:
    if not (cocos_index in range(1, 9) or cocos_index in range(11, 19)):
        raise ValueError(
            f"Convention number {cocos_index} is not valid. "
            "Must be between 1 and 8 or 11 and 18."
        )
    # gets first matching index
    return next(x for x in all_conventions() if x.cc_index == cocos_index)


def all_conventions() -> tuple[COCOSConvention, ...]:
    return (
        C1,
        C2,
        C3,
        C4,
        C5,
        C6,
        C7,
        C8,
        C11,
        C12,
        C13,
        C14,
        C15,
        C16,
        C17,
        C18,
    )


def matching_conventions(
    exp_Bp: Optional[ZeroOne] = None,
    sign_Bp: Optional[Sign] = None,
    sign_R_phi_Z: Optional[Sign] = None,
    sign_rho_theta_phi: Optional[Sign] = None,
    sign_q: Optional[Sign] = None,
    sign_pprime: Optional[Sign] = None,
) -> list[COCOSConvention]:
    def _match_cocos(c: COCOSConvention) -> bool:
        return all(
            [
                exp_Bp is None or c.exp_Bp == exp_Bp,
                sign_Bp is None or c.sign_Bp == sign_Bp,
                sign_R_phi_Z is None or c.sign_R_phi_Z == sign_R_phi_Z,
                sign_rho_theta_phi is None
                or c.sign_rho_theta_phi == sign_rho_theta_phi,
            ]
        )

    return list(filter(_match_cocos, all_conventions()))


def identify_eqdsk(
    eqdsk: EQDSKInterface,
    # is the direction of the toroidal B field
    # (the toroidal magnetic field) clockwise?
    clockwise_phi: Optional[bool] = None,
    volt_seconds_per_radian: Optional[bool] = None,
    flux_surfaces_minor_radius: Optional[np.ndarray] = None,
) -> list[COCOSConvention]:
    sign_Ip = np.sign(eqdsk.cplasma)
    sign_B0 = np.sign(eqdsk.bcentre)

    # todo: what if qpsi is None??
    sign_q = np.sign(eqdsk.qpsi)
    if sign_q.min() != sign_q.max():
        raise ValueError("The sign of qpsi is not consistent across the flux surfaces.")
    sign_q = sign_q.max()

    sign_d_psi_towards_boundary = np.sign(eqdsk.psibdry - eqdsk.psimag)

    sign_Bp = sign_Ip * sign_d_psi_towards_boundary
    sign_rho_theta_phi = sign_Ip * sign_B0 * sign_q

    if clockwise_phi is None:
        sign_R_phi_Z = None
    elif clockwise_phi:
        sign_R_phi_Z = Sign.N
    else:
        sign_R_phi_Z = Sign.P

    # not sure if this is needed
    # in the utf8 code they do it, in omas they don't
    # if sign_R_phi_Z is not None:
    #     sign_rho_theta_phi *= sign_R_phi_Z.value

    sign_Bp = Sign(sign_Bp)
    sign_rho_theta_phi = Sign(sign_rho_theta_phi)

    exp_Bp = None
    if volt_seconds_per_radian is not None:
        if volt_seconds_per_radian:
            exp_Bp = ZeroOne.ONE
        else:
            exp_Bp = ZeroOne.ZERO
    elif flux_surfaces_minor_radius is not None:
        a = flux_surfaces_minor_radius

        # from OMAS
        # https://gafusion.github.io/omas/_modules/omas/omas_physics.html

        # idx of min q, not sure why yet
        index = np.argmin(np.abs(eqdsk.qpsi))
        if index == 0:
            index = 1

        q_estimate = np.abs(
            (np.pi * eqdsk.bcentre * (a[index] - a[0]) ** 2)
            / (
                # todo: this won't work because psi is a 2d grid?
                eqdsk.psi[index]
                - eqdsk.psi[0]
            )
        )
        q_actual = np.abs(eqdsk.qpsi[index])

        if abs(q_estimate - q_actual) < abs(q_estimate / (2 * np.pi) - q_actual):
            exp_Bp = ZeroOne.ONE
        else:
            exp_Bp = ZeroOne.ZERO

    return matching_conventions(
        exp_Bp=exp_Bp,
        sign_Bp=sign_Bp,
        sign_R_phi_Z=sign_R_phi_Z,
        sign_rho_theta_phi=sign_rho_theta_phi,
    )


def convert_eqdsk(eqdsk: EQDSKInterface, to_cocos_index: int) -> EQDSKInterface:
    org_cocos = eqdsk.cocos_convention

    if org_cocos.cc_index == to_cocos_index:
        return eqdsk

    tgt_cocos = get_convention(to_cocos_index)

    org_eqdsk = eqdsk
    tgt_eqdsk = deepcopy(eqdsk)

    eff_bp = tgt_cocos.sign_Bp.value * org_cocos.sign_Bp.value
    eff_R_phi_Z = tgt_cocos.sign_R_phi_Z.value * org_cocos.sign_R_phi_Z.value
    eff_rho_theta_phi = (
        tgt_cocos.sign_rho_theta_phi.value * org_cocos.sign_rho_theta_phi.value
    )
    eff_exp_bp = tgt_cocos.exp_Bp.value - org_cocos.exp_Bp.value

    # when eff_exp_bp is -1, it means org is vs/rad and trg isn't,
    # thus this is 1/2pi.
    # Meaning for tgt_psi, we'll effectively multiply org_psi by 2pi
    # getting rid of the /rad factor.
    # pprime and ffprime go with 1/psi thus are the opposite.
    pi_factor = 2 * np.pi**eff_exp_bp

    tgt_eqdsk.cplasma = eff_R_phi_Z * org_eqdsk.cplasma
    tgt_eqdsk.bcentre = eff_R_phi_Z * org_eqdsk.bcentre

    # todo: not sure about these, but this make sense to me
    # tgt_eqdsk.Ic = eff_R_phi_Z * org_eqdsk.Ic
    # tgt_eqdsk.fpol = eff_R_phi_Z * eff_rho_theta_phi * org_eqdsk.fpol

    tgt_eqdsk.psi = eff_bp * eff_R_phi_Z * (1 / pi_factor) * org_eqdsk.psi
    tgt_eqdsk.psibdry = eff_bp * eff_R_phi_Z * (1 / pi_factor) * org_eqdsk.psibdry
    tgt_eqdsk.psimag = eff_bp * eff_R_phi_Z * (1 / pi_factor) * org_eqdsk.psimag

    tgt_eqdsk.pprime = eff_bp * eff_R_phi_Z * pi_factor * org_eqdsk.pprime
    tgt_eqdsk.ffprime = eff_bp * eff_R_phi_Z * pi_factor * org_eqdsk.ffprime

    if org_eqdsk.qpsi is not None:
        # there isn't much agreement on this one
        tgt_eqdsk.qpsi = eff_R_phi_Z * eff_rho_theta_phi * org_eqdsk.qpsi

    tgt_eqdsk.identify()

    return tgt_eqdsk
