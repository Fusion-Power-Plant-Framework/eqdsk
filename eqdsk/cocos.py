from dataclasses import dataclass
from enum import Enum
from typing import Optional

import numpy as np

# from eqdsk.error import COCOSError
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
    sign_q: Sign  # same as sign_rho_theta_phi
    sign_pprime: Sign  # opposite of sign_Bp


C1 = COCOSConvention(
    cc_index=1,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.N,
)
C11 = COCOSConvention(
    cc_index=11,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.N,
)
C2 = COCOSConvention(
    cc_index=2,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.N,
)
C12 = COCOSConvention(
    cc_index=12,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.N,
)
C3 = COCOSConvention(
    cc_index=3,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.P,
)
C13 = COCOSConvention(
    cc_index=13,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.P,
)
C4 = COCOSConvention(
    cc_index=4,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.P,
)
C14 = COCOSConvention(
    cc_index=14,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.P,
)
C5 = COCOSConvention(
    cc_index=5,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.N,
)
C15 = COCOSConvention(
    cc_index=15,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.N,
)
C6 = COCOSConvention(
    cc_index=6,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.N,
)
C16 = COCOSConvention(
    cc_index=16,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.P,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.N,
    sign_q=Sign.N,
    sign_pprime=Sign.N,
)
C7 = COCOSConvention(
    cc_index=7,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.P,
)
C17 = COCOSConvention(
    cc_index=17,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.P,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.P,
)
C8 = COCOSConvention(
    cc_index=8,
    exp_Bp=ZeroOne.ZERO,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.P,
)
C18 = COCOSConvention(
    cc_index=18,
    exp_Bp=ZeroOne.ONE,
    sign_Bp=Sign.N,
    sign_R_phi_Z=Sign.N,
    sign_rho_theta_phi=Sign.P,
    sign_q=Sign.P,
    sign_pprime=Sign.P,
)


class COCOSHelper:
    @staticmethod
    def get_convention(cocos_index: int) -> COCOSConvention:
        if not (cocos_index in range(1, 9) or cocos_index in range(11, 19)):
            raise ValueError(
                f"Convention number {cocos_index} is not valid. "
                "Must be between 1 and 8 or 11 and 18."
            )
        return next(
            x for x in COCOSHelper.all_conventions() if x.cc_index == cocos_index
        )

    @staticmethod
    def all_conventions() -> list[COCOSConvention]:
        return [
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
        ]

    @staticmethod
    def matching(
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
                    sign_q is None or c.sign_q == sign_q,
                    sign_pprime is None or c.sign_pprime == sign_pprime,
                ]
            )

        return list(filter(_match_cocos, COCOSHelper.all_conventions()))

    @staticmethod
    def identify_eqdsk(
        eqdsk: EQDSKInterface,
        # is the direction of the toroidal B field
        # (the toroidal magnetic field) clockwise?
        clockwise_phi: Optional[bool] = None,
        flux_surfaces_minor_radius: Optional[np.ndarray] = None,
    ) -> list[COCOSConvention]:
        sign_Ip = np.sign(eqdsk.cplasma)
        sign_B0 = np.sign(eqdsk.bcentre)

        sign_q = np.sign(eqdsk.qpsi)
        if sign_q.min() != sign_q.max():
            raise ValueError(
                "The sign of qpsi is not consistent across the flux surfaces."
            )
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
        if flux_surfaces_minor_radius is not None:
            a = flux_surfaces_minor_radius
            # from OMAS
            # https://gafusion.github.io/omas/_modules/omas/omas_physics.html
            index = np.argmin(np.abs(eqdsk.qpsi))
            if index == 0:
                index = 1
            q_estimate = np.abs(
                (np.pi * eqdsk.bcentre * (a[index] - a[0]) ** 2)
                / (
                    # this won't work becuase psi is a 2d grid
                    eqdsk.psi[index]
                    - eqdsk.psi[0]
                )
            )
            q_actual = np.abs(eqdsk.qpsi[index])
            if abs(q_estimate - q_actual) < abs(q_estimate / (2 * np.pi) - q_actual):
                exp_Bp = ZeroOne.ONE
            else:
                exp_Bp = ZeroOne.ZERO

        return COCOSHelper.matching(
            exp_Bp=exp_Bp,
            sign_Bp=sign_Bp,
            sign_R_phi_Z=sign_R_phi_Z,
            sign_rho_theta_phi=sign_rho_theta_phi,
        )

    @staticmethod
    def convert_eqdsk(eqdsk: EQDSKInterface, to_cocos_index: int) -> EQDSKInterface:
        ...


# def convert(
#     eqdsk: EQDSKInterface, conv_to: COCOS, conv_from: Optional[COCOS] = None
# ) -> EQDSKInterface:
#     conv_from = conv_from or identify(eqdsk)
#     if isinstance(conv_from, Iterable):
#         if len(conv_from) == 1:
#             conv_from = conv_from[0]
#         else:
#             raise COCOSError("More than one COCOS version available: {conv_from}")
#     if isinstance(conv_from, COCOS) and conv_from == conv_to:
#         eqdsk_warn("No conversion needed")
#         return eqdsk

#     eqdsk_print("Converting from COCOS-{conv_from.number} to COCOS-{conv_to.number}")

#     new_eqdsk = deepcopy(eqdsk)
#     update_dict = {}

#     # I have no idea what I'm doing for the rest of this function...
#     two_pi_exp = (np.pi * 2) ** (conv_to.exp_Bp.value - conv_from.exp_Bp.value)
#     sign_Ip = conv_from.sign_R_phi_Z.value * conv_to.sign_R_phi_Z.value
#     sign_Bp = conv_from.sign_Bp.value * conv_to.sign_Bp.value
#     sign_rtp = conv_from.sign_rho_theta_phi.value * conv_to.sign_rho_theta_phi.value

#     if int(sign_Ip * sign_Bp * two_pi_exp) != 1:
#         transform = int(sign_Ip * sign_Bp) * two_pi_exp
#         transform2 = int(sign_Ip * sign_Bp) / two_pi_exp
#         update_dict["psi"] = new_eqdsk.psi * transform
#         update_dict["psibdry"] = new_eqdsk.psibdry * transform
#         update_dict["psimag"] = new_eqdsk.psimag * transform
#         update_dict["pprime"] = new_eqdsk.pprime * transform2
#         update_dict["ffprime"] = new_eqdsk.ffprime * transform2

#     if int(sign_Ip * sign_Ip * sign_rtp) != 1:
#         update_dict["qpsi"] = new_eqdsk.qpsi * int(sign_Ip * sign_Ip * sign_rtp)

#     if int(sign_Ip) != 1:
#         update_dict["bcentre"] = new_eqdsk.bcentre * sign_Ip
#         update_dict["Ic"] = new_eqdsk.Ic * sign_Ip

#     if int(sign_Ip * sign_rtp) != 1:
#         update_dict["fpol"] = new_eqdsk.fpol * sign_Ip * sign_rtp
#         # Missing one...
#         # transforms['BP'] = sign_Ip * sign_rtp

#     new_eqdsk.update(update_dict)

#     return new_eqdsk
