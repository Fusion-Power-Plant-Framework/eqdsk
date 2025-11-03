# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""Read and write EQDSK files."""

from __future__ import annotations

import json
import os
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import fortranformat as ff
import numpy as np

from eqdsk.cocos import COCOS, KnownCOCOS, convert_eqdsk, identify_eqdsk
from eqdsk.errors import (
    MissingQpsiDataError,
    NoSingleConventionError,
)
from eqdsk.log import eqdsk_print, eqdsk_warn
from eqdsk.models import Sign
from eqdsk.tools import is_num, json_writer

if TYPE_CHECKING:
    from collections.abc import Iterator, Sized
    from io import TextIOWrapper

    import numpy.typing as npt


EQDSK_EXTENSIONS = [".eqdsk", ".eqdsk_out", ".geqdsk"]


@dataclass
class EQDSKInterface:
    """Container for data from an EQDSK file.

    The G-EQDSK file format is described here:
        https://fusion.gat.com/conferences/snowmass/working/mfe/physics/p3/equilibria/g_eqdsk_s.pdf

    Notes
    -----
    G-EQDSK is from the 1980's and EQDSK files should generally only be
    read and not written. New equilibria should really just be saved as
    JSON files.
    """

    DEFAULT_COCOS: ClassVar[int] = 11

    bcentre: float
    """Vacuum toroidal Magnetic field at the reference radius [T]."""
    cplasma: float
    """Plasma current [A]."""
    dxc: np.ndarray
    """X half-thicknesses of the coils [m]."""
    dzc: np.ndarray
    """Z half-thicknesses of the coils [m]."""
    ffprime: np.ndarray
    """FF' function on 1-D flux grid [m.T^2/V.s/rad]."""
    fpol: np.ndarray
    """Poloidal current function f = R*B on 1-D flux [T.m]."""
    Ic: np.ndarray
    """Coil currents [A]."""
    name: str
    """Name of the equilibrium EQDSK [dimensionless]."""
    nbdry: int
    """Number of boundary points [dimensionless]."""
    ncoil: int
    """Number of coils [dimensionless]."""
    nlim: int
    """Number of limiters [dimensionless]."""
    nx: int
    """Number of grid points in the radial direction [dimensionless]."""
    nz: int
    """Number of grid points in the vertical direction [dimensionless]."""
    pprime: np.ndarray
    """P' function on 1-D flux grid [N/m^2/V.s/rad]."""
    pressure: np.ndarray
    """Plasma pressure function on 1-D flux grid [N/m^2]."""
    psi: np.ndarray
    """Poloidal magnetic flux on the 2-D grid
    [V.s/rad (COCOS 1-8) or V.s (COCOS 11-18)]."""
    psibdry: float
    """Poloidal flux at the plasma boundary (LCFS) [V.s/rad]."""
    psimag: float
    """Poloidal flux at the magnetic axis [V.s/rad]."""
    xbdry: np.ndarray
    """X coordinates of the plasma boundary [m]."""
    xc: np.ndarray
    """X coordinates of the coils [m]."""
    xcentre: float
    """Radius of the reference toroidal magnetic [m]."""
    xdim: float
    """Horizontal dimension of the spatial grid [m]."""
    xgrid1: float
    """Minimum radius of the spatial grid [m]."""
    xlim: np.ndarray
    """X coordinates of the limiters [m]."""
    xmag: float
    """Radius of the magnetic axis [m]."""
    zbdry: np.ndarray
    """Z coordinates of the plasma boundary [m]."""
    zc: np.ndarray
    """Z coordinates of the coils [m]."""
    zdim: float
    """Vertical dimension of the spatial grid [m]."""
    zlim: np.ndarray
    """Z coordinates of the limiters [m]."""
    zmag: float
    """Z coordinate of the magnetic axis [m]."""
    zmid: float
    """Z coordinate of the middle of the spatial grid [m]."""
    x: np.ndarray | None = None
    """X 1-D vector [m] (calculated if not given)."""
    z: np.ndarray | None = None
    """Z 1-D vector [m] (calculated if not given)."""
    psinorm: np.ndarray | None = None
    """Normalised psi vector [A] (calculated if not given)."""
    qpsi: np.ndarray | None = None
    """Safety factor values on the 1-D flux grid [dimensionless]."""
    file_name: str | None = None
    """The EQDSK file the data originates from."""
    coil_names: list[str] | None = None
    """Name of the coils"""
    coil_types: list[str] | None = None
    """Type of the coils"""
    comment: str | None = None
    """Any comment stored on file"""
    unprocessed_data: np.ndarray | None = None
    """Any unprocessed data from raw eqdsks"""

    def __post_init__(self):
        """Calculate derived parameters if they're not given."""
        if self.x is None:
            self.x = _derive_x(self.xgrid1, self.xdim, self.nx)
        if self.z is None:
            self.z = _derive_z(self.zmid, self.zdim, self.nz)
        if self.psinorm is None:
            self.psinorm = _derive_psinorm(self.fpol)
        self._cocos = None

    def __repr__(self) -> str:  # noqa: D105
        if self.qpsi is None:
            conv_p = ", ".join(
                str(c.index) for c in identify_eqdsk(self, qpsi_positive=True)
            )
            conv_n = ", ".join(
                str(c.index) for c in identify_eqdsk(self, qpsi_positive=False)
            )
            conventions_str = f"{conv_p} (+ve qpsi) | {conv_n} (-ve qpsi)"
        else:
            conventions_str = ", ".join(str(c.index) for c in identify_eqdsk(self))

        if self._cocos is not None:
            conventions_str = f"{conventions_str} (current COCOS: {self._cocos.index})"

        return f"""
Identifiable COCOS: {conventions_str}

Main properties:

* bcentre: {self.bcentre} [T] (vacuum toroidal magnetic field at the reference radius)
* cplasma: {self.cplasma} [A] (plasma current)

* xcentre: {self.xcentre} [m] (radius of the reference toroidal magnetic axis)

* xmag: {self.xmag} [m] (radius of the magnetic axis)
* zmag: {self.zmag} [m] (z coordinate of the magnetic axis)

* psibdry: {self.psibdry} [V.s/rad or V.s] (poloidal flux at the plasma boundary (LCFS))
* psimag: {self.psimag} [V.s/rad or V.s] (poloidal flux at the magnetic axis)

* nlim: {self.nlim} (number of limiters)

File properties:

* name: {self.name}
* file_name: {self.file_name}

Grid properties:

* nx: {self.nx} (no. points in the radial direction)
* nz: {self.nz} (no. points in the vertical direction)

* xdim: {self.xdim} [m] (grid radial extent)
* zdim: {self.zdim} [m] (grid vertical extent)
"""

    @classmethod
    def from_file(
        cls,
        file_path: str | Path,
        from_cocos: int | str | COCOS | KnownCOCOS | None = None,
        to_cocos: int | str | COCOS | KnownCOCOS | None = DEFAULT_COCOS,
        *,
        clockwise_phi: bool | None = None,
        volt_seconds_per_radian: bool | None = None,
        qpsi_positive: bool | None = None,
        no_cocos: bool = False,
    ) -> EQDSKInterface:
        """Create an EQDSKInterface object from a file.

        Parameters
        ----------
        file_path:
            Path to a file of one of the following formats:
                - JSON
                - eqdsk
                - eqdsk_out
                - geqdsk
        from_cocos:
            The user set COCOS of the EQDSK file.
            This sets what COCCOS the file will be id's as
            and will raise it's not one of the determined COCOS.
        to_cocos:
            The COCOS to convert the EQDSK file to. If None, the file
            will not be converted.
        clockwise_phi:
            Whether the EQDSK file's phi is clockwise or not.
        volt_seconds_per_radian:
            Whether the EQDSK file's psi is in volt seconds per radian.
        no_cocos:
            Whether to return the EQDSK data without performing
            any identifying COCOS identification or conversion.
        qpsi_positive:
            Whether qpsi is positive or not, required for identification
            when qpsi is not present in the file.

        Returns
        -------
        :
            An instance of this class containing the EQDSK file's data.

        Raises
        ------
        ValueError
            Unknown file format
        NoSingleConventionError
            More than one COCOS convention identified
        """
        file_path = Path(file_path)
        file_extension = file_path.suffix
        file_name = file_path.name
        if file_extension.lower() in EQDSK_EXTENSIONS:
            inst = cls(file_name=file_name, **_read_eqdsk(file_path))
        elif file_extension.lower() == ".json":
            inst = cls(file_name=file_name, **_read_json(file_path))
        else:
            raise ValueError(f"Unrecognised file format '{file_extension}'.")

        if no_cocos:
            return inst

        try:
            inst.identify(
                as_cocos=from_cocos,
                clockwise_phi=clockwise_phi,
                volt_seconds_per_radian=volt_seconds_per_radian,
                qpsi_positive=qpsi_positive,
            )
        except NoSingleConventionError as e:
            raise NoSingleConventionError(
                e.conventions,
                message_extra="You need to specify `from_cocos` or "
                "`clockwise_phi` and `volt_seconds_per_radian`.",
            ) from None

        if to_cocos is not None:
            inst = inst.to_cocos(to_cocos)

        return inst

    @property
    def cocos(self) -> COCOS:
        """Return the COCOS for this eqdsk.

        Raises
        ------
        ValueError
            COCOS not identified
        """
        if self._cocos is None:
            raise ValueError(
                "The COCOS for this eqdsk has not yet been identified. "
                "The 'identify' method must be called first "
                "before the COCOS can be returned.",
            )
        return self._cocos

    def identify(
        self,
        as_cocos: int | str | COCOS | KnownCOCOS | None = None,
        *,
        clockwise_phi: bool | None = None,
        volt_seconds_per_radian: bool | None = None,
        qpsi_positive: bool | None = None,
    ):
        """Identifies the COCOS of this eqdsk.

        Note
        ----
        This sets the internal _cocos attribute and does not return
        anything.

        Parameters
        ----------
        as_cocos:
            The COCOS index to convert the EQDSK file to.
            If given, the file will be id's as the given COCOS,
            only if it one of the possible identified COCOS.
        clockwise_phi:
            Whether the EQDSK file's phi is clockwise or not.
        volt_seconds_per_radian:
            Whether the EQDSK file's psi is in volt seconds per radian.
        qpsi_positive:
            Whether qpsi is positive or not, required for identification
            when self.qpsi is None.

        Raises
        ------
        ValueError
            If as_cocos is given but does not match any identified COCOS.
        ValueError
            If no COCOS can be identified.
        MissingQpsiDataError
            qpsi not provided or found in file and qpsi_positive is not set.
        """
        qpsi_sign = None if qpsi_positive is None else Sign(qpsi_positive)
        qpsi_is_not_set = self.qpsi is None or np.allclose(self.qpsi, 0)

        if qpsi_is_not_set:
            if qpsi_sign:
                eqdsk_warn(
                    "eqdsk contains no qpsi data, but "
                    f"`qpsi_positive={qpsi_positive}` provided. "
                    f"Setting qpsi to array of {qpsi_sign.value}'s."
                )
                self.qpsi = np.ones(self.nx) * qpsi_sign.value
            else:
                raise MissingQpsiDataError(
                    message_extra="To resolve this, set the `qpsi_positive` parameter. "
                    "This is the sign (true: 1, false: -1) "
                    "of qpsi across the flux surfaces.\n"
                    "If you are unsure what the sign should be, "
                    "refer to the COCOS spec (or its implementation "
                    "in this package) and see which standard fits the direction"
                    "for theta and phi (CW or CCW) for this EQDSK.\n"
                    "You can also experiment by setting `qpsi_positive` and checking if "
                    "the resulting COCOS('s) is(are) correct.",
                )

        conventions = identify_eqdsk(
            self,
            clockwise_phi=clockwise_phi,
            volt_seconds_per_radian=volt_seconds_per_radian,
        )

        def _id():
            if as_cocos:
                cocos_fmt = COCOS(as_cocos)
                matching_conv = next((c for c in conventions if c == cocos_fmt), None)
                if not matching_conv:
                    raise ValueError(
                        f"No convention found that matches "
                        f"the given COCOS index {cocos_fmt.index}, "
                        f"from the possible ({', '.join([str(c.index) for c in conventions])}).",  # noqa: E501
                    )
                return matching_conv
            if len(conventions) != 1:
                raise NoSingleConventionError(
                    conventions,
                    message_extra="You need to specify `as_cocos` or "
                    "`clockwise_phi` and `volt_seconds_per_radian`.",
                )
            return conventions[0]

        c = _id()
        eqdsk_print(f"EQDSK identified as COCOS {c.index}.")
        self._cocos = c

    def to_cocos(self, to_cocos: int | str | COCOS | KnownCOCOS) -> EQDSKInterface:
        """
        Returns
        -------
        :
             A copy of this eqdsk converted to the given COCOS.

        Note
        ----
        This returns a new instance of the EQDSKInterface class.
        """
        to_cocos = COCOS(to_cocos)
        if self.cocos == to_cocos:
            return self
        eqdsk_print(f"Converting EQDSK to COCOS {to_cocos.index}.")
        return convert_eqdsk(self, to_cocos.index)

    def to_dict(self, *, with_comment: bool = False) -> dict:
        """
        Returns
        -------
        :
            A dictionary of the EQDSK data.
        """
        d = asdict(self)
        # Remove the file name as this is metadata, not EQDSK data
        del d["file_name"]
        del d["unprocessed_data"]
        if not with_comment:
            d.pop("comment")

        return d

    def write(
        self,
        file_path: str | Path,
        file_format: str = "json",
        json_kwargs: dict | None = None,
        *,
        strict_spec: bool = True,
        write_comment: bool = False,
    ):
        """Write the EQDSK data to file in the given format.

        Parameters
        ----------
        file_path:
            Path to where the file should be written.
        file_format:
            The format to save the file in. One of 'json', 'eqdsk', or
            'geqdsk'.
        json_kwargs:
            Key word arguments to pass to the ``json.dump`` call. Only
            used if ``format`` is 'json'.
        strict_spec:
            As https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf arrays have the format
            5e16.9, disabling this changes the format to 5ES23.16e2
        write_comment:
            write any comments to file
        """
        if file_format == "json":
            json_kwargs = {} if json_kwargs is None else json_kwargs
            json_writer(
                self.to_dict(with_comment=write_comment), file_path, **json_kwargs
            )
        elif file_format in {"eqdsk", "geqdsk"}:
            eqdsk_warn(
                "You are in the 21st century. "
                "Are you sure you want to be making an EDQSK in this day and age?"
            )
            _write_eqdsk(
                file_path,
                self.to_dict(with_comment=write_comment),
                strict_spec=strict_spec,
            )

    def update(self, eqdsk_data: dict[str, Any]):
        """Update this object's data with values from a dictionary.

        Parameters
        ----------
        eqdsk_data:
            A dict containing the new eqdsk data.

        Raises
        ------
        ValueError
            If a key in `eqdsk_data` does not correspond to an
            attribute of this class.
        """
        for key, value in eqdsk_data.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise ValueError(
                    f"Cannot update EQDSKInterface from dict. Unrecognised key '{key}'.",
                )


def _read_json(file_path: Path) -> dict[str, Any]:
    with file_path.open() as file:
        data = json.load(file)

    for k, value in data.items():
        if isinstance(value, list) and k not in {"coil_type", "coil_names"}:
            data[k] = np.asarray(value)

    # For backward compatibility where 'psinorm' was sometimes 'pnorm'
    if "pnorm" in data:
        if "psinorm" in data:
            del data["pnorm"]
        else:
            data["psinorm"] = data.pop("pnorm")

    return data


def _read_array(
    tokens: Iterator[str], n: int, name: str = "Unknown"
) -> npt.NDArray[float]:
    data = np.zeros([n])
    try:
        for i in np.arange(n):
            data[i] = float(next(tokens))
    except StopIteration as e:
        raise ValueError(f"Failed reading array {name} of size {n}") from e
    return data


def _read_2d_array(
    tokens: Iterator[str], n_x: int, n_y: int, name: str = "Unknown"
) -> npt.NDArray[float]:
    data = np.zeros([n_y, n_x])
    for i in np.arange(n_y):
        data[i, :] = _read_array(tokens, n_x, f"{name}[{i!s}]]")
    return np.transpose(data)


def _eqdsk_generator(file: TextIOWrapper) -> Iterator[str]:
    """Transform a file object into a generator, following G-EQDSK number
    conventions.

    Parameters
    ----------
    file:
        The file to read

    Yields
    ------
    :
        The lines of the file
    """
    while line := file.readline():
        yield from _process_exponents(line).split()


def _process_exponents(line: str) -> str:
    """Distinguish negative/positive numbers from negative/positive exponent"""  # noqa: DOC201
    if "E" in line or "e" in line:
        sym = "__SYM__"  # Symbol to avoid replacing valid text
        line = line.replace("E-", sym)
        line = line.replace("e-", sym)
        line = line.replace("-", " -")
        line = line.replace(sym, "e-")
        line = line.replace("E+", sym)
        line = line.replace("e+", sym)
        line = line.replace("+", " ")
        line = line.replace(sym, "e+")
    return line


CHAR_IN_LINE = re.compile(r"^[0-9 .e\-+]+$")


def _get_line(file: TextIOWrapper) -> str:
    position = file.tell()
    while line := file.readline():
        if line.strip():
            break
    file.seek(position)
    return _process_exponents(line)


def _coils_check(file: TextIOWrapper) -> bool:
    line = _get_line(file)
    sline = line.split()
    if sline and CHAR_IN_LINE.match(line):
        try:
            int(sline[0])
        except ValueError:
            return False
        else:
            return True
    return False


def _more_data_check(file: TextIOWrapper) -> bool:
    line = _get_line(file)
    return bool(line.split() and CHAR_IN_LINE.match(line))


def _eqdsk_out_of_spec_generator(
    file: TextIOWrapper, *, more_data: bool = False
) -> Iterator[str]:
    if more_data:
        while line := file.readline():
            yield from _process_exponents(line).split()
    else:
        while line := file.readline():
            yield line


def _read_eqdsk(file_path: Path) -> dict:
    with file_path.open("r") as file:
        description = file.readline()
        if not description:
            raise OSError(f"Could not read the file '{file}'.")
        description = description.split()

        ints = [v for v in description if is_num(v)]
        if len(ints) < 3:  # noqa: PLR2004
            raise OSError(
                f"Should be at least 3 numbers in the first line of the EQDSK {file}.",
            )

        data = {}
        n_x = int(ints[-2])
        n_z = int(ints[-1])
        data["name"] = description[0]
        data["nx"] = n_x
        data["nz"] = n_z

        tokens = _eqdsk_generator(file)
        for name in [
            "xdim",
            "zdim",
            "xcentre",
            "xgrid1",
            "zmid",
            "xmag",
            "zmag",
            "psimag",
            "psibdry",
            "bcentre",
            "cplasma",
            "psimag",
            None,
            "xmag",
            None,
            "zmag",
            None,
            "psibdry",
            None,
            None,
        ]:
            if name is not None:  # Lots of dummies and duplication
                data[name] = float(next(tokens))
            else:
                next(tokens)  # Dummy

        for name in ["fpol", "pressure", "ffprime", "pprime"]:
            data[name] = _read_array(tokens, n_x, name)

        data["psi"] = _read_2d_array(tokens, n_x, n_z, "psi")
        data["qpsi"] = _read_array(tokens, n_x, "qpsi")
        if np.allclose(data["qpsi"], 0):
            data["qpsi"] = None
        nbdry = int(next(tokens))
        nlim = int(next(tokens))
        data["nbdry"] = nbdry
        data["nlim"] = nlim

        xbdry = np.zeros(nbdry)
        zbdry = np.zeros(nbdry)
        for i in range(nbdry):
            xbdry[i] = float(next(tokens))
            zbdry[i] = float(next(tokens))
        data["xbdry"] = xbdry
        data["zbdry"] = zbdry

        xlim = np.zeros(nlim)
        zlim = np.zeros(nlim)
        for i in range(nlim):
            xlim[i] = float(next(tokens))
            zlim[i] = float(next(tokens))
        data["xlim"] = xlim
        data["zlim"] = zlim

        # Additional utility data
        data["x"] = _derive_x(data["xgrid1"], data["xdim"], data["nx"])
        data["z"] = _derive_z(data["zmid"], data["zdim"], data["nz"])
        data["psinorm"] = _derive_psinorm(data["fpol"])

        # end of spec everything past this point is to support custom extensions
        has_coils = _coils_check(file)
        oos_tokens = _eqdsk_out_of_spec_generator(file, more_data=has_coils)

        data["ncoil"] = int(next(oos_tokens)) if has_coils else 0
        data["xc"], data["zc"], data["dxc"], data["dzc"], data["Ic"] = (
            _get_coils_from_eqdsk(data["ncoil"], oos_tokens)
        )

        data["unprocessed_data"] = (
            _get_extra_data(oos_tokens) if _more_data_check(file) else None
        )

        comments = _get_comment(_eqdsk_out_of_spec_generator(file, more_data=False))
        if comments:
            data["comment"] = comments

    return data


def _get_coils_from_eqdsk(ncoil: int, tokens: Iterator[str]) -> tuple[list[float], ...]:
    x_c = np.zeros(ncoil)
    z_c = np.zeros(ncoil)
    dxc = np.zeros(ncoil)
    dzc = np.zeros(ncoil)
    i_c = np.zeros(ncoil)
    for i in range(ncoil):
        x_c[i] = float(next(tokens))
        z_c[i] = float(next(tokens))
        dxc[i] = float(next(tokens))
        dzc[i] = float(next(tokens))
        i_c[i] = float(next(tokens))
    return x_c, z_c, dxc, dzc, i_c


def _get_extra_data(tokens: Iterator[str]) -> npt.NDArray[float]:
    data = []
    try:
        while d := float(next(tokens)):
            data.append(d)
    except (ValueError, StopIteration):
        pass
    return np.array(data)


def _get_comment(tokens: Iterator[str]) -> str:
    comments = []
    try:
        while True:
            comments.append(next(tokens))
    except StopIteration:
        pass
    return "".join(comments).strip("\n")


def _derive_x(xgrid1: float, xdim: float, nx: int) -> npt.NDArray[float]:
    return np.linspace(xgrid1, xgrid1 + xdim, nx)


def _derive_z(zmid: float, zdim: float, nz: int) -> npt.NDArray[float]:
    return np.linspace(zmid - zdim / 2, zmid + zdim / 2, nz)


def _derive_psinorm(fpol: Sized) -> npt.NDArray[float]:
    return np.linspace(0, 1, len(fpol))


def _write_eqdsk(file_path: str | Path, data: dict, *, strict_spec: bool = True):  # noqa: PLR0915
    """Write data out to a text file in G-EQDSK format.

    Parameters
    ----------
    file_path:
        The full path string of the file to be created
    data:
        Dictionary of EQDSK data.
    strict_spec:
        As https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf arrays have the format
        5e16.9, disabling this changes the format to 5ES23.16e2

    Raises
    ------
    ValueError
        Length of qpsi is not consistent with the number of grid points
    """
    file_path = Path(file_path)
    if file_path.suffix not in EQDSK_EXTENSIONS:
        file_path = file_path.with_suffix(".eqdsk")

    with Path(file_path).open("w") as file:

        def write_header(
            fortran_format: ff.FortranRecordWriter,
            id_string: str,
            var_list: list[str],
        ):
            """Write a G-EQDSK header out to file.

            Parameters
            ----------
            fortran_format:
                FortranRecordWriter object for Fortran format edit
                descriptor to be used for header output.
            id_string:
                String containing name of file to be used as identification
                string. Will be trimmed if length exceeds 39 characters,
                so it will fit within the permitted header length of the
                GEQDSK specification when a timestamp is added.
            var_list:
                List of names of keys in EQDSKInterface.data identifying
                variables to add to the header following the id_string.
                Empty strings will be recorded as 0.
            """
            line = [id_string]
            line += [data[v] if v else 0 for v in var_list]
            file.write(fortran_format.write(line))
            file.write("\n")

        def write_line(fortran_format: ff.FortranRecordWriter, var_list: list[str]):
            """Write a line of variable values out to a G-EQDSK file.

            Parameters
            ----------
            fortran_format:
                FortranRecordWriter object for Fortran format edit
                descriptor to be used for the format of the line output.
            var_list:
                List of names of keys in EQDSKInterface.data identifying
                variables to added to the current line.
                Empty strings will be recorded as 0.
            """
            line = [data[v] if v else 0 for v in var_list]
            file.write(fortran_format.write(line))
            file.write("\n")

        def write_array(fortran_format: ff.FortranRecordWriter, array: np.ndarray):
            """Write a numpy array out to a G-EQDSK file.

            Parameters
            ----------
            fortran_format:
                FortranRecordWriter object for Fortran format edit
                descriptor to be used for the format of the line output.
            array:
                Numpy array of variables to be written to file.
                Array will be flattened in column-major (Fortran)
                order if is more than one-dimensional.
            """
            if array.ndim > 1:
                flat_array = array.flatten(order="F")
                file.write(fortran_format.write(flat_array))
            else:
                file.write(fortran_format.write(array))
            file.write("\n")

        # Create id string for file comprising of timestamp and trimmed filename
        # that fits the 48 character limit of strings in EQDSK headers.
        timestamp = time.strftime("%d%m%Y")
        trimmed_name = data["name"][0 : 48 - len(timestamp) - 1]
        file_id_string = f"{trimmed_name}_{timestamp}"

        # Define dummy data for qpsi if it has not been previously defined.
        qpsi = (
            np.zeros(data["nx"]) if data["qpsi"] is None else np.atleast_1d(data["qpsi"])
        )
        if np.allclose(qpsi, 0):
            eqdsk_warn(
                f"Writing EQDSK with qpsi all zeros. {data['name']}, {data['qpsi']}"
            )
        if len(qpsi) == 1:
            qpsi = np.full(data["nx"], qpsi)
        elif len(qpsi) != data["nx"]:
            raise ValueError(
                "the length of qpsi should be 1 or the number of x grid points"
            )

        # Create array containing coilset information.
        coil = np.zeros(5 * data["ncoil"])
        for i, value in enumerate(["xc", "zc", "dxc", "dzc", "Ic"]):
            coil[i::5] = data[value]

        # Create FortranRecordWriter objects with the Fortran format
        # edit descriptors to be used in the G-EQDSK output.
        f2000 = ff.FortranRecordWriter("a48,3i4")
        f2020 = ff.FortranRecordWriter("5e16.9" if strict_spec else "5ES23.16e2")
        f2022 = ff.FortranRecordWriter("2i5")
        fCSTM = ff.FortranRecordWriter("i5")

        # Write header in f2000 (6a8,3i4) format.
        write_header(f2000, file_id_string, ["", "nx", "nz"])
        # Write out lines containing floats in f2020 (5e16.9) format.
        write_line(f2020, ["xdim", "zdim", "xcentre", "xgrid1", "zmid"])
        write_line(f2020, ["xmag", "zmag", "psimag", "psibdry", "bcentre"])
        write_line(f2020, ["cplasma", "psimag", "", "xmag", ""])
        write_line(f2020, ["zmag", "", "psibdry", "", ""])
        # Write out arrays in in f2020 (5e16.9) format.
        write_array(f2020, data["fpol"])
        write_array(f2020, data["pressure"])
        write_array(f2020, data["ffprime"])
        write_array(f2020, data["pprime"])
        write_array(f2020, data["psi"])
        write_array(f2020, qpsi)
        # Write out number of boundary points and limiters f2022 (2i5) format.
        write_line(f2022, ["nbdry", "nlim"])
        # Write out boundary point and limiter data as array of ordered pairs.
        write_array(f2020, np.array([data["xbdry"], data["zbdry"]]))
        write_array(f2020, np.array([data["xlim"], data["zlim"]]))

        # Output of coilset information. This is an extension to the
        # regular eqdsk format.
        if data["ncoil"] > 0:
            write_line(fCSTM, ["ncoil"])
            write_array(f2020, coil)

        if comment := data.get("comment"):
            cl = list(filter(lambda x: x.strip(), comment.split("\n")))[1:]
            if len(cl) > 1:
                comment_char = os.path.commonprefix(cl) or " " * 4
                if not comment.startswith(comment_char):
                    comment = f"{comment_char}{comment}"
            file.write(f"\n{comment}\n")
