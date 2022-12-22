from copy import deepcopy
from dataclasses import dataclass
from enum import Enum, EnumMeta
from types import DynamicClassAttribute
from typing import Iterable, Optional, Union

import numpy as np

from eqdsk.error import COCOSError
from eqdsk.file import EQDSKInterface
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

    exp_Bp: expBp
    sign_Bp: Bp
    sign_R_phi_Z: RpZ
    sign_rho_theta_phi: rtp
    sign_q: q
    sign_pprime: pp


class COCOSEnumMeta(EnumMeta):
    """
    Allow override of KeyError error string
    """

    def __getitem__(self, name: Union[float, int, str]) -> Enum:
        if isinstance(name, (float, int)):
            return super().__getitem__(f"C{int(name)}")
        else:
            try:
                return super().__getitem__(name)
            except KeyError:
                return super().__getitem__(f"C{name}")


class COCOS(Enum, metaclass=COCOSEnumMeta):
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

    @DynamicClassAttribute
    def number(self) -> int:
        return int(self.name.strip("C"))

    @DynamicClassAttribute
    def exp_Bp(self) -> int:
        return self.value.exp_Bp

    @DynamicClassAttribute
    def sign_R_phi_Z(self) -> int:
        return self.value.sign_R_phi_Z

    @DynamicClassAttribute
    def sign_rho_theta_phi(self) -> int:
        return self.value.sign_rho_theta_phi

    @DynamicClassAttribute
    def sign_q(self) -> int:
        return self.value.sign_q

    @DynamicClassAttribute
    def sign_pprime(self) -> int:
        return self.value.sign_pprime


def convert(
    eqdsk: EQDSKInterface, conv_to: COCOS, conv_from: Optional[COCOS] = None
) -> EQDSKInterface:
    conv_from = conv_from or identify(eqdsk)
    if isinstance(conv_from, Iterable):
        if len(conv_from) == 1:
            conv_from = conv_from[0]
        else:
            raise COCOSError("More than one COCOS version available: {conv_from}")
    if isinstance(conv_from, COCOS) and conv_from == conv_to:
        eqdsk_warn("No conversion needed")
        return eqdsk

    eqdsk_print("Converting from COCOS-{conv_from.number} to COCOS-{conv_to.number}")

    new_eqdsk = deepcopy(eqdsk)
    update_dict = {}

    # I have no idea what I'm doing here...
    two_pi_exp = np.pi * 2**conv_to.exp_Bp - conv_from.exp_Bp
    sign_Ip = conv_from.sign_R_phi_Z * conv_to.sign_R_phi_Z
    sign_Bp = conv_from.sign_Bp * conv_to.sign_Bp
    sign_rtp = conv_from.sign_rho_theta_phi * conv_to.sign_rho_theta_phi

    if int(sign_Ip * sign_Bp * two_pi_exp) != 1:
        transform = int(sign_Ip * sign_Bp) * two_pi_exp
        transform2 = int(sign_Ip * sign_Bp) / two_pi_exp
        update_dict["psi"] = new_eqdsk.psi * transform
        update_dict["psibdry"] = new_eqdsk.psibdry * transform
        update_dict["psimag"] = new_eqdsk.psimag * transform
        update_dict["pprime"] = new_eqdsk.pprime * transform2
        update_dict["ffprime"] = new_eqdsk.ffprime * transform2

    if int(sign_Ip * sign_Ip * sign_rtp) != 1:
        update_dict["qpsi"] = new_eqdsk.qpsi * int(sign_Ip * sign_Ip * sign_rtp)

    if int(sign_Ip) != 1:

        update_dict["bcentre"] = new_eqdsk.bcentre * sign_Ip
        update_dict["Ic"] = new_eqdsk.Ic * sign_Ip

    if int(sign_Ip * sign_rtp) != 1:
        update_dict["fpol"] = new_eqdsk.fpol * sign_Ip * sign_rtp
        # Missing one...
        # transforms['BP'] = sign_Ip * sign_rtp

    eqdsk.update(update_dict)


def identify(
    eqdsk: EQDSKInterface,
    sign_R_phi_Z: Optional[int] = None,
    minor_radius_fl: Optional[np.ndarray] = None,
) -> COCOS:
    current_sign = np.sign(eqdsk.Ic)

    if sign_R_phi_Z is None:
        sign_R_phi_Z = [RpZ(1), RpZ(-1)]
    else:
        sign_R_phi_Z = [RpZ(int(sign_R_phi_Z))]

    raw_number = []
    for s_RpZ in sign_R_phi_Z:

        sign_Bz = -1 * s_RpZ.value * current_sign
        sign_Bp = -1 * sign_Bz * s_RpZ.value * np.sign(eqdsk.psibdry - eqdsk.psimag)
        sign_rho_theta_phi = (
            s_RpZ.value * current_sign * np.sign(eqdsk.bcentre) * np.sign(eqdsk.qpsi)
        )
        raw_number.append(
            COCOS._identification[Bp(int(sign_Bp))][rtp(int(sign_rho_theta_phi))][s_RpZ]
        )

    raw_number = np.array(raw_number)

    if minor_radius_fl:
        index = np.argmin(np.abs(eqdsk.qpsi)) or 1
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

    return (COCOS[no] for no in number)


@dataclass(repr=False)
class COCOSEQDSKInterface(EQDSKInterface):
    @property
    def cocos_version(self):
        # TODO set sign_R_phi_Z and minor_radius_fl or have this as a raw attr
        return identify(self)

    def as_cocos(self, cocos_conv: Union[float, int, str]):
        return convert(self, COCOS[cocos_conv])  # TODO save cocos version on object
