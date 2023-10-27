from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from enum import Enum, unique
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from eqdsk.file import EQDSKInterface


class ZeroOne(Enum):
    ZERO = 0
    ONE = 1


class Sign(Enum):
    POSITIVE = 1
    NEGATIVE = -1


@dataclass(frozen=True)
class COCOSValues:
    cc_index: int
    exp_Bp: ZeroOne  # 0 -> ψ/2pi^(1-0) i.e. ψ/2pi, 1 -> ψ/2pi^(1-1) i.e. just ψ
    sign_Bp: Sign  # ψ (P -> increasing)/(N -> decreasing) from the magnetic axis
    sign_R_phi_Z: Sign  # P -> (R, ϕ, Z) ϕ cnt-clockwise, N -> (R, Z, ϕ) ϕ clockwise
    sign_rho_theta_phi: Sign  # P -> (ρ, θ, ϕ), N -> (ρ, ϕ, θ)


@unique
class COCOS(Enum):
    C1 = COCOSValues(
        cc_index=1,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C11 = COCOSValues(
        cc_index=11,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C2 = COCOSValues(
        cc_index=2,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C12 = COCOSValues(
        cc_index=12,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C3 = COCOSValues(
        cc_index=3,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C13 = COCOSValues(
        cc_index=13,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C4 = COCOSValues(
        cc_index=4,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C14 = COCOSValues(
        cc_index=14,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C5 = COCOSValues(
        cc_index=5,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C15 = COCOSValues(
        cc_index=15,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C6 = COCOSValues(
        cc_index=6,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C16 = COCOSValues(
        cc_index=16,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.POSITIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.NEGATIVE,
    )
    C7 = COCOSValues(
        cc_index=7,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C17 = COCOSValues(
        cc_index=17,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.POSITIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C8 = COCOSValues(
        cc_index=8,
        exp_Bp=ZeroOne.ZERO,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )
    C18 = COCOSValues(
        cc_index=18,
        exp_Bp=ZeroOne.ONE,
        sign_Bp=Sign.NEGATIVE,
        sign_R_phi_Z=Sign.NEGATIVE,
        sign_rho_theta_phi=Sign.POSITIVE,
    )

    def __init__(self, c: COCOSValues):
        # shortcuts .value access
        self.cc_index = c.cc_index
        self.exp_Bp = c.exp_Bp
        self.sign_Bp = c.sign_Bp
        self.sign_R_phi_Z = c.sign_R_phi_Z
        self.sign_rho_theta_phi = c.sign_rho_theta_phi

    @classmethod
    def with_index(cls, cocos_index: int) -> COCOSValues:
        if not (cocos_index in range(1, 9) or cocos_index in range(11, 19)):
            raise ValueError(
                f"Convention number {cocos_index} is not valid. "
                "Must be between 1 and 8 or 11 and 18."
            )
        return next(x for x in cls if x.cc_index == cocos_index)

    @classmethod
    def matching_convention(
        cls,
        exp_Bp: ZeroOne,
        sign_Bp: Sign,
        sign_R_phi_Z: Sign,
        sign_rho_theta_phi: Sign,
    ) -> COCOSValues:
        def _match_cocos(c: COCOS) -> bool:
            return all(
                [
                    c.exp_Bp == exp_Bp,
                    c.sign_Bp == sign_Bp,
                    c.sign_R_phi_Z == sign_R_phi_Z,
                    c.sign_rho_theta_phi == sign_rho_theta_phi,
                ]
            )

        return next(filter(_match_cocos, cls))


def identify_eqdsk(
    eqdsk: EQDSKInterface,
    clockwise_phi: Optional[bool] = None,
    volt_seconds_per_radian: Optional[bool] = None,
) -> list[COCOS]:
    if eqdsk.qpsi is None:
        raise ValueError("qpsi is not defined in the eqdsk file.")

    conventions = []
    for cw_phi in [True, False] if clockwise_phi is None else [clockwise_phi]:
        for vs_pr in (
            [True, False]
            if volt_seconds_per_radian is None
            else [volt_seconds_per_radian]
        ):
            conventions.append(
                identify_cocos(
                    plasma_current=eqdsk.cplasma,
                    b_center=eqdsk.bcentre,
                    psi_at_boundary=eqdsk.psibdry,
                    psi_at_mag_axis=eqdsk.psimag,
                    q_psi=eqdsk.qpsi,
                    phi_clockwise_from_top=cw_phi,
                    volt_seconds_per_radian=vs_pr,
                )
            )
    # return sort asc by cc_index
    conventions.sort(key=lambda x: x.cc_index)
    return conventions


def identify_cocos(
    plasma_current: float,
    b_center: float,
    psi_at_boundary: float,
    psi_at_mag_axis: float,
    q_psi: np.ndarray,
    phi_clockwise_from_top: bool,
    volt_seconds_per_radian: bool,
) -> COCOS:
    if phi_clockwise_from_top:
        sign_R_phi_Z = Sign.NEGATIVE
    else:
        sign_R_phi_Z = Sign.POSITIVE

    if volt_seconds_per_radian:
        exp_Bp = ZeroOne.ONE
    else:
        exp_Bp = ZeroOne.ZERO

    sign_Ip = np.sign(plasma_current)
    sign_B0 = np.sign(b_center)

    sign_q = np.sign(q_psi)
    if sign_q.min() != sign_q.max():
        raise ValueError("The sign of qpsi is not consistent across the flux surfaces.")
    sign_q = sign_q.max()

    sign_d_psi_towards_boundary = np.sign(psi_at_boundary - psi_at_mag_axis)

    sign_Bp = sign_Ip * sign_d_psi_towards_boundary
    sign_rho_theta_phi = sign_Ip * sign_B0 * sign_q * sign_R_phi_Z.value

    sign_Bp = Sign(sign_Bp)
    sign_rho_theta_phi = Sign(sign_rho_theta_phi)

    return COCOS.matching_convention(
        exp_Bp=exp_Bp,
        sign_Bp=sign_Bp,
        sign_R_phi_Z=sign_R_phi_Z,
        sign_rho_theta_phi=sign_rho_theta_phi,
    )


def convert_eqdsk(eqdsk: EQDSKInterface, to_cocos_index: int) -> EQDSKInterface:
    org_cocos = eqdsk.cocos
    tgt_cocos = COCOS.with_index(to_cocos_index)

    org_eqdsk = eqdsk
    tgt_eqdsk = deepcopy(eqdsk)

    if org_cocos.cc_index == tgt_cocos.cc_index:
        return tgt_eqdsk

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
    tgt_eqdsk.Ic = eff_R_phi_Z * org_eqdsk.Ic
    tgt_eqdsk.fpol = eff_R_phi_Z * eff_rho_theta_phi * org_eqdsk.fpol

    tgt_eqdsk.psi = eff_bp * eff_R_phi_Z * (1 / pi_factor) * org_eqdsk.psi
    tgt_eqdsk.psibdry = eff_bp * eff_R_phi_Z * (1 / pi_factor) * org_eqdsk.psibdry
    tgt_eqdsk.psimag = eff_bp * eff_R_phi_Z * (1 / pi_factor) * org_eqdsk.psimag

    tgt_eqdsk.pprime = eff_bp * eff_R_phi_Z * pi_factor * org_eqdsk.pprime
    tgt_eqdsk.ffprime = eff_bp * eff_R_phi_Z * pi_factor * org_eqdsk.ffprime

    if org_eqdsk.qpsi is not None:
        tgt_eqdsk.qpsi = eff_R_phi_Z * eff_rho_theta_phi * org_eqdsk.qpsi

    # not sure if this should be done here, but kinda makes sense
    # because you have the tgt_cocos
    tgt_eqdsk.identify(
        clockwise_phi=tgt_cocos.sign_R_phi_Z == Sign.NEGATIVE,
        volt_seconds_per_radian=tgt_cocos.exp_Bp == ZeroOne.ONE,
    )

    return tgt_eqdsk
