from dataclasses import dataclass
from enum import Enum

import numpy as np

from eqdsk.log import eqdsk_print, eqdsk_warn


class expBp(Enum):

    ZERO = 0
    ONE = 1


class Bp(Enum):

    P = 1
    N = -1


class RpZ(Enum):

    P = 1
    N = -1


class rtp(Enum):

    P = 1
    N = -1


class q(Enum):

    P = 1
    N = -1


class pp(Enum):

    P = 1
    N = -1


@dataclass
class Convention:

    exp_Bp: exp_Bp
    sign_Bp: Bp
    sign_R_phi_Z: RpZ
    sign_rho_theta_phi: rtp
    sign_q: q
    sign_pprime: pp


class COCOS(Enum):
    C1: Convention = Convention(expBp(0), Bp(1), RpZ(1), rtp(1), q(1), pp(-1))
    C2: Convention = Convention(expBp(0), Bp(1), RpZ(-1), rtp(1), q(1), pp(-1))
    C3: Convention = Convention(expBp(0), Bp(-1), RpZ(1), rtp(-1), q(-1), pp(1))
    C4: Convention = Convention(expBp(0), Bp(-1), RpZ(-1), rtp(-1), q(-1), pp(1))
    C5: Convention = Convention(expBp(0), Bp(1), RpZ(1), rtp(-1), q(-1), pp(-1))
    C6: Convention = Convention(expBp(0), Bp(1), RpZ(-1), rtp(-1), q(-1), pp(-1))
    C7: Convention = Convention(expBp(0), Bp(-1), RpZ(1), rtp(1), q(1), pp(1))
    C8: Convention = Convention(expBp(0), Bp(-1), RpZ(-1), rtp(1), q(1), pp(1))
    C11: Convention = Convention(expBp(1), Bp(1), RpZ(1), rtp(1), q(1), pp(-1))
    C12: Convention = Convention(expBp(1), Bp(1), RpZ(-1), rtp(1), q(1), pp(-1))
    C13: Convention = Convention(expBp(1), Bp(-1), RpZ(1), rtp(-1), q(-1), pp(1))
    C14: Convention = Convention(expBp(1), Bp(-1), RpZ(-1), rtp(-1), q(-1), pp(1))
    C15: Convention = Convention(expBp(1), Bp(1), RpZ(1), rtp(-1), q(-1), pp(-1))
    C16: Convention = Convention(expBp(1), Bp(1), RpZ(-1), rtp(-1), q(-1), pp(-1))
    C17: Convention = Convention(expBp(1), Bp(-1), RpZ(1), rtp(1), q(1), pp(1))
    C18: Convention = Convention(expBp(1), Bp(-1), RpZ(-1), rtp(1), q(1), pp(1))

    _identification = {
        Bp(1): {rtp(1): {RpZ(1): 1, RpZ(-1): 2}, rtp(-1): {RpZ(1): 5, RpZ(-1): 6}},
        Bp(-1): {rtp(-1): {RpZ(1): 3, RpZ(-1): 4}, rtp(1): {RpZ(1): 7, RpZ(-1): 8}},
    }

    def number(self):
        return self.name.strip("C")


def convert(
    eqdsk: EQDSKInterface, conv_to: COCOS, conv_from: Optional[COCOS] = None
) -> EQDSKInterface:
    conv_from = conv_from or identify(eqdsk)
    if isinstance(conv_from, Iterable) and len(conv_from) == 1:
        conv_from = conv_from[0]
    if isinstance(conv_from, COCOS) and conv_from == conv_to:
        eqdsk_warn("No conversion needed")
        return eqdsk

    eqdsk_print("Converting from COCOS-{conv_from.number} to COCOS-{conv_to.number}")

    # eqdsk.


def identify(eqdsk: EQDSKInterface, sign_R_phi_Z=None, minor_radius_fl=None) -> COCOS:
    current_sign = np.sign(eqdsk.Ic)

    if sign_R_phi_Z is None:
        sign_R_phi_Z = [1, -1]
    else:
        sign_R_phi_Z = [sign_R_phi_Z]

    for s_R_phi_Z in sign_R_phi_Z:

        sign_Bz = -1 * s_R_phi_Z * current_sign
        sign_Bp = sign_Bz * s_R_phi_Z * -1 * np.sign(eqdsk.psibdry - eqdsk.psimag)
        sign_rho_theta_phi = (
            s_R_phi_Z * current_sign * np.sign(eqdsk.bcentre) * np.sign(eqdsk.qpsi)
        )
        raw_number.append(
            COCOS._identification[Bp(sign_Bp)][rtp(sign_rho_theta_phi)][RpZ(s_R_phi_Z)]
        )

    raw_number = np.array(raw_number)

    if minor_radius_fl:
        index = np.argmin(np.abs(eqdsk.qpsi))
        if index == 0:
            index = 1
        q_estimate = np.abs(
            (np.pi * eqdsk.bcentre * (minor_radius_fl - minor_radius_fl[0]) ** 2)
            / (eqdsk.psimag - eqdsk.psimag[0])
        )
        number = raw_number + 10 * int(
            np.abs(q_estimate[index] - eqdsk.qpsi[index])
            < np.abs(q_estimate[index] / (2 * np.pi) - eqdsk.qpsi[index])
        )
    else:
        number = np.append(raw_number, raw_number + 10)

    return (COCOS[f"C{no}"] for no in number)
