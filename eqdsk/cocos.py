from dataclasses import dataclass

import numpy as np

from eqdsk.log import eqdsk_warn, eqdsk_print


@dataclass
class Convention:

    number: int
    exp_Bp: int
    sign_Bp: int
    sign_R_phi_Z: int
    sign_rho_theta_phi: int
    sign_q: int
    sign_pprime: int

@dataclass
class COCOS:
    C1 : Convention = Convention(1, 0, 1, 1, 1, 1, -1)
    C2 : Convention = Convention(2, 0, 1, -1, 1, 1, -1)
    C3 : Convention = Convention(3, 0, -1, 1, -1, -1, 1)
    C4 : Convention = Convention(4, 0, -1, -1, -1, -1, 1)
    C5 : Convention = Convention(5, 0, 1, 1, -1, -1, -1)
    C6 : Convention = Convention(6, 0, 1, -1, -1, -1, -1)
    C7 : Convention = Convention(7, 0, -1, 1, 1, 1, 1)
    C8 : Convention = Convention(8, 0, -1, -1, 1, 1, 1)
    C11 : Convention = Convention(11, 1, 1, 1, 1, 1, -1)
    C12 : Convention = Convention(12, 1, 1, -1, 1, 1, -1)
    C13 : Convention = Convention(13, 1, -1, 1, -1, -1, 1)
    C14 : Convention = Convention(14, 1, -1, -1, -1, -1, 1)
    C15 : Convention = Convention(15, 1, 1, 1, -1, -1, -1)
    C16 : Convention = Convention(16, 1, 1, -1, -1, -1, -1)
    C17 : Convention = Convention(17, 1, -1, 1, 1, 1, 1)
    C18 : Convention = Convention(18, 1, -1, -1, 1, 1, 1)

    _identification = {1:{-1:{1:5,-1:6},1:{1:1,-1:2}} , -1:{-1:{1:3,-1:4},1:{1:7,-1:8}}}


def convert(eqdsk: EQDSKInterface, conv_to: Convention, conv_from: Optional[Convention] = None) -> EQDSKInterface:
    conv_from = conv_from or identify(eqdsk)
    if isinstace(conv_from, Convention) or len(conv_from) == 1:
        conv_from = conv_from[0]
        if conv_from == conv_to:
            eqdsk_warn("No conversion needed")
            return eqdsk

    eqdsk_print("Converting from COCOS-{conv_from.number} to COCOS-{conv_to.number}")

    eqdsk.


def identify(eqdsk: EQDSKInterface, sign_R_phi_Z=None, minor_radius_fl=None) -> Convention:
    current_sign = np.sign(eqdsk.Ic)

    if sign_R_phi_Z is None:
        sign_R_phi_Z = [1, -1]
    else:
        sign_R_phi_Z = [sign_R_phi_Z]

    for s_R_phi_Z in sign_R_phi_Z:

        sign_Bz = -1*s_R_phi_Z* current_sign
        sign_Bp = sign_Bz* s_R_phi_Z* -1 * np.sign(eqdsk.psibdry - eqdsk.psimag)
        sign_rho_theta_phi = s_R_phi_Z * current_sign *np.sign(eqdsk.bcentre) * np.sign(eqdsk.qpsi)
        raw_number.append(COCOS._identification[sign_Bp][sign_rho_theta_phi][s_R_phi_Z])

    raw_number = np.array(raw_number)

    if minor_radius_fl:
        index = np.argmin(np.abs(eqdsk.qpsi))
        if index == 0:
            index = 1
        q_estimate = np.abs((np.pi * eqdsk.bcentre * (minor_radius_fl - minor_radius_fl[0]) ** 2) / (eqdsk.psimag - eqdsk.psimag[0]))
        number = raw_number + 10*int(np.abs(q_estimate[index] - eqdsk.qpsi[index]) < np.abs(q_estimate[index] / (2 * np.pi) - eqdsk.qpsi[index]))
    else
        number = np.append(raw_number, raw_number+10)

    return getattr(COCOS, f"C{no}") for no in number

