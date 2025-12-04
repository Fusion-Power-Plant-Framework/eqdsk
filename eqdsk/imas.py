# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""IMAS integration for EQDSK."""

from __future__ import annotations

from contextlib import suppress
from typing import TYPE_CHECKING

import imas
import numpy as np
from imas.exception import DataEntryException, IDSNameError

if TYPE_CHECKING:
    from imas.ids_primitive import IDSPrimitive

    from eqdsk.file import EQDSKInterface


def _unwrap_imas_value(value: IDSPrimitive, /, *, default=None):
    if value.has_value:
        return value.value

    return default


def from_imas(
    db: imas.DBEntry,
    time_index: int = 0,
    profiles_2d_index: int = 0,
    time: float | None = None,
) -> dict:
    try:
        equilibrium_top_level = db.get("equilibrium")
    except IDSNameError as e:
        raise RuntimeError(
            "'equilibrium' top-level IDS is not present in the database"
        ) from e

    # It should not be an issue if the limiter is not in the database
    # handle it by not returning limiter quantities
    try:
        limiter = db.get("wall").description_2d[0].limiter.unit[0].outline

        xlim = _unwrap_imas_value(limiter.r)
        zlim = _unwrap_imas_value(limiter.z)
    except (DataEntryException, IDSNameError):
        xlim = np.array([])
        zlim = np.array([])

    time_index = (
        time_index
        if time is None
        else np.argmin(np.abs(equilibrium_top_level.time.value - time))
    )
    eq_time = equilibrium_top_level.time_slice[time_index]
    eq_vtf = equilibrium_top_level.vacuum_toroidal_field
    global_quantities = eq_time.global_quantities
    boundary_outline = eq_time.boundary.outline
    eq_profiles_1d = eq_time.profiles_1d
    eq_profiles_2d = eq_time.profiles_2d[profiles_2d_index]

    gridx = eq_profiles_2d.grid.dim1.value
    gridz = eq_profiles_2d.grid.dim2.value
    xdim = max(gridx) - min(gridx)
    zdim = max(gridz) - min(gridz)

    if hasattr(eq_vtf, "b0"):
        xcentre = eq_vtf.r0.value
        bcentre = eq_vtf.b0[time_index]
    else:
        xcentre = (max(xdim) + min(xdim)) / 2
        bcentre = (
            global_quantities.magnetic_axis.b_field_phi.value
            * global_quantities.magnetic_axis.r.valie
            / xcentre
        )

    dxc = []
    dzc = []
    xc = []
    zc = []
    ic = []
    coil_names = []
    coil_types = []

    # Again, its not critical if we cannot access the PF coil
    # data so provide sensible empty entries.
    with suppress(DataEntryException, IDSNameError):
        pf_top_level = db.get("pf_active")

        for coil in pf_top_level.coil:
            if coil.geometry.geometry_type.value == 2:
                dxc.append(coil.geometry.rectangle.width.value / 2)
                dzc.append(coil.geometry.rectangle.height.value / 2)
                xc.append(coil.geometry.rectangle.width.value)
                zc.append(coil.geometry.rectangle.height.value)
                ic.append(coil.current.data)
                coil_names.append(coil.name.value)
                coil_types.append("PF")

    ncoil = len(coil_names)

    return {
        "bcentre": bcentre,
        "cplasma": _unwrap_imas_value(global_quantities.ip),
        "dxc": np.array(dxc),
        "dzc": np.array(dzc),
        "ffprime": _unwrap_imas_value(eq_profiles_1d.f_df_dpsi),
        "fpol": _unwrap_imas_value(eq_profiles_1d.f, default=np.array([])),
        "Ic": np.array(ic),
        "name": _unwrap_imas_value(
            equilibrium_top_level.ids_properties.name, default=""
        ),
        "file_name": db.uri,
        "nbdry": len(boundary_outline.r),
        "ncoil": ncoil,
        "nlim": len(xlim),
        "nx": gridx.size,
        "nz": gridz.size,
        "pprime": _unwrap_imas_value(eq_profiles_1d.dpressure_dpsi),
        "pressure": _unwrap_imas_value(eq_profiles_1d.pressure),
        "psi": _unwrap_imas_value(eq_profiles_2d.psi),
        "psibdry": _unwrap_imas_value(global_quantities.psi_boundary),
        "psimag": _unwrap_imas_value(
            global_quantities.psi_axis,
            default=_unwrap_imas_value(global_quantities.psi_magnetic_axis),
        ),
        "xbdry": _unwrap_imas_value(boundary_outline.r),
        "xc": np.array(xc),
        "xcentre": xcentre,
        "xdim": xdim,
        "xgrid1": min(gridx),
        "xlim": xlim,
        "xmag": _unwrap_imas_value(global_quantities.magnetic_axis.r),
        "zbdry": _unwrap_imas_value(boundary_outline.z),
        "zc": np.array(zc),
        "zdim": zdim,
        "zlim": zlim,
        "zmag": _unwrap_imas_value(global_quantities.magnetic_axis.z),
        "zmid": (max(gridz) + min(gridz)) / 2,
        "qpsi": _unwrap_imas_value(eq_profiles_1d.q),
        "coil_names": coil_names,
        "coil_types": coil_types,
    }


# TODO!
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
