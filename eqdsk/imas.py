# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""IMAS integration for EQDSK."""

from __future__ import annotations

from contextlib import contextmanager
from enum import Enum
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from eqdsk.file import EQDSKInterface

IMAS = MappingProxyType({
    "cplasma": "global_quantities.ip",
    # "dxc":
    # "dzc":
    "ffprime": "profiles_1d.f_df_dpsi",
    "fpol": "profiles_1d.f",
    # "Ic":
    "pprime": "profiles_1d.dpressure_dpsi",
    "pressure": "profiles_1d.pressure",
    "psi": "psi",
    "psibdry": "global_quantities.psi_boundary",
    "psimag": "global_quantities.psi_axis",
    "xbdry": "boundary.outline.r",
    # "xc":
    "xlim": "wall.description_2d.0.limiter.unit.0.outline.r",
    "xmag": "global_quantities.magnetic_axis.r",
    "zbdry": "boundary.outline.z",
    # "zc":
    "zlim": "wall.description_2d.0.limiter.unit.0.outline.z",
    "zmag": "global_quantities.magnetic_axis.z",
    "qpsi": "profiles_1d.q",
    # "coil_names":
    # "coil_types":
})


class ReadWrite(Enum):
    READ = "r"
    WRITE = "w"


@contextmanager
def imas_connection(path: str, mode: ReadWrite = ReadWrite.READ) -> DBEntry:
    import imaspy

    with imaspy.DBEntry(f"imas:hdf5?path={path}", mode.value()) as db:
        yield db


def from_imas(
    db, time_index: int = 0, profiles_2d_index: int = 0, time: float | None = None
) -> dict:
    time_index = (
        time_index if time is None else np.argmin(np.abs(db["equilibrium.time"] - time))
    )
    eq_time = db[f"equilibrium.time_slice.{time_index}"]
    eq_vtf = db["equilibrium.vacuum_toroidal_field"]
    eq_profiles_2d = eq_time[f"profiles_2d.{profiles_2d_index}"]
    gridx = eq_profiles_2d["grid.dim1"]
    gridz = eq_profiles_2d["grid.dim2"]
    xdim = max(gridx) - min(gridx)
    zdim = max(gridz) - min(gridz)
    xlim = db[IMAS["xlim"]]
    xbdry = db[IMAS["xbdry"]]
    if "b0" in eq_vtf:
        xcentre = eq_vtf["r0"]
        bcentre = eq_vtf["b0"][time_index]
    else:
        xcentre = (max(xdim) + min(xdim)) / 2
        bcentre = (
            db["global_quantities.magnetic_axis.b_field_tor"]
            * db[IMAS["xmag"]]
            / xcentre
        )

    return {
        "bcentre": bcentre,
        "cplasma": db[IMAS["cplasma"]],
        # "dxc": ,
        # "dzc": ,
        "ffprime": db[IMAS["ffprime"]],
        "fpol": db[IMAS["fpol"]],
        # "Ic": ,
        "name": " TODO shot info here",
        "nbdry": len(xbdry),
        # "ncoil": ,
        "nlim": len(xlim),
        "nx": gridx.size,
        "nz": gridz.size,
        "pprime": db[IMAS["pprime"]],
        "pressure": db[IMAS["pressure"]],
        "psi": eq_profiles_2d[IMAS["psi"]],
        "psibdry": db[IMAS["psibdry"]],
        "psimag": db[IMAS["psimag"]],
        "xbdry": xbdry,
        # "xc": ,
        "xcentre": xcentre,
        "xdim": xdim,
        "xgrid1": min(gridx),
        "xlim": xlim,
        "xmag": db[IMAS["xmag"]],
        "zbdry": db[IMAS["xbdry"]],
        # "zc": ,
        "zdim": zdim,
        "zlim": db[IMAS["zlim"]],
        "zmag": db[IMAS["zmag"]],
        "zmid": (max(gridz) + min(gridz)) / 2,
        "qpsi": db[IMAS["qpsi"]],
        # "coil_names": ,
        # "coil_types": ,
    }


def to_imas(db, eqdsk: EQDSKInterface, time_index: int = 0):
    eqdsk = eqdsk.to_cocos(1)

    db["dataset_description.data_entry.pulse"] = 0
    db["equilibrium.ids_properties.comment"] = "eqdsk python package"
    db[IMAS["cplasma"]] = eqdsk.cplasma
    db[IMAS["ffprime"]] = eqdsk.ffprime
    db[IMAS["fpol"]] = eqdsk.fpol
    db[IMAS["pprime"]] = eqdsk.pprime
    db[IMAS["pressure"]] = eqdsk.pressure
    db[IMAS["psibdry"]] = eqdsk.psibdry
    db[IMAS["psimag"]] = eqdsk.psimag
    db[IMAS["xmag"]] = eqdsk.xmag
    db[IMAS["xbdry"]] = eqdsk.xbdry
    db[IMAS["zlim"]] = eqdsk.zlim
    db[IMAS["zmag"]] = eqdsk.zmag
    db[IMAS["qpsi"]] = eqdsk.qpsi

    eq_time = db[f"equilibrium.time_slice.{time_index}"]
    eq_time["time"] = 0
    eq_profiles_2d = eq_time["profiles_2d.0"]
    eq_profiles_2d[IMAS["psi"]] = eqdsk.psi
    eq_profiles_2d["grid_type.index"] = 1
    eq_profiles_2d["grid.dim1"] = np.linspace(0, eqdsk.xdim, eqdsk.nx) + eqdsk.xgrid1
    eq_profiles_2d["grid.dim2"] = (
        np.linspace(0, eqdsk.zdim, eqdsk.nz) - eqdsk.zdim / 2
    ) + eqdsk.zmid

    eq_vtf = db["equilibrium.vacuum_toroidal_field"]
    eq_vtf["r0"] = eqdsk.xcentre
    eq_vtf["b0"][time_index] = eqdsk.bcentre

    wdb = db["wall.description_2d.0.limiter"]
    wdb["type.name"] = "first_wall"
    wdb["type.index"] = 0
    wdb["type.description"] = "first wall"
    wdb["unit.0.outline.r"] = eqdsk.xlim
    wdb["unit.0.outline.z"] = eqdsk.zlim
    db["wall.time"][time_index] = eq_time["time"]
