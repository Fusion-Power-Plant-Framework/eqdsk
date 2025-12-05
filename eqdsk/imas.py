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

        xlim = _unwrap_imas_value(limiter.r, default=np.array([]))
        zlim = _unwrap_imas_value(limiter.z, default=np.array([]))
    except (DataEntryException, IDSNameError):
        xlim = np.array([])
        zlim = np.array([])

    time_index = (
        time_index
        if time is None
        else np.argmin(np.abs(_unwrap_imas_value(equilibrium_top_level.time) - time))
    )
    eq_time = equilibrium_top_level.time_slice[time_index]
    eq_vtf = equilibrium_top_level.vacuum_toroidal_field
    global_quantities = eq_time.global_quantities
    boundary_outline = eq_time.boundary.outline
    eq_profiles_1d = eq_time.profiles_1d
    eq_profiles_2d = eq_time.profiles_2d[profiles_2d_index]

    gridx = _unwrap_imas_value(eq_profiles_2d.grid.dim1)
    gridz = _unwrap_imas_value(eq_profiles_2d.grid.dim2)
    xdim = max(gridx) - min(gridx)
    zdim = max(gridz) - min(gridz)

    if hasattr(eq_vtf, "b0"):
        xcentre = _unwrap_imas_value(eq_vtf.r0)
        bcentre = eq_vtf.b0[time_index]
    else:
        xcentre = (max(xdim) + min(xdim)) / 2
        bcentre = (
            _unwrap_imas_value(global_quantities.magnetic_axis.b_field_phi)
            * _unwrap_imas_value(global_quantities.magnetic_axis.r)
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
            if _unwrap_imas_value(coil.geometry.geometry_type) == 2:
                dxc.append(_unwrap_imas_value(coil.geometry.rectangle.width) / 2)
                dzc.append(_unwrap_imas_value(coil.geometry.rectangle.height) / 2)
                xc.append(_unwrap_imas_value(coil.geometry.rectangle.width))
                zc.append(_unwrap_imas_value(coil.geometry.rectangle.height))
                ic.append(coil.current.data)
                coil_names.append(_unwrap_imas_value(coil.name))
                coil_types.append("PF")

    ncoil = len(coil_names)

    psibdry = _unwrap_imas_value(global_quantities.psi_boundary)
    psimag = _unwrap_imas_value(
        global_quantities.psi_axis,
        default=_unwrap_imas_value(global_quantities.psi_magnetic_axis),
    )
    psinorm = _unwrap_imas_value(eq_profiles_1d.psi_norm)
    if psinorm is None and psibdry is not None and psinorm is not None:
        psi1d = _unwrap_imas_value(eq_profiles_1d.psi)
        if psi1d is not None:
            psinorm = (psi1d - psimag) / (psibdry - psimag)

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
        "psinorm": psinorm,
        "psi": _unwrap_imas_value(eq_profiles_2d.psi),
        "psibdry": psibdry,
        "psimag": psimag,
        "xbdry": _unwrap_imas_value(boundary_outline.r, default=np.array([])),
        "xc": np.array(xc),
        "xcentre": xcentre,
        "xdim": xdim,
        "xgrid1": min(gridx),
        "xlim": xlim,
        "xmag": _unwrap_imas_value(global_quantities.magnetic_axis.r),
        "zbdry": _unwrap_imas_value(boundary_outline.z, default=np.array([])),
        "zc": np.array(zc),
        "zdim": zdim,
        "zlim": zlim,
        "zmag": _unwrap_imas_value(global_quantities.magnetic_axis.z),
        "zmid": (max(gridz) + min(gridz)) / 2,
        "qpsi": _unwrap_imas_value(eq_profiles_1d.q),
        "coil_names": coil_names,
        "coil_types": coil_types,
    }


def to_imas(
    db: imas.DBEntry,
    eqdsk: EQDSKInterface,
    time_index: int = 0,
    profiles_2d_index: int = 0,
    ids_factory_kwargs=None,
):
    eqdsk = eqdsk.to_cocos(1)

    ids_factory = imas.IDSFactory(**(ids_factory_kwargs or {}))

    equilibrium_ids = ids_factory.equilibrium()
    vacuum_toroidal_field = equilibrium_ids.vacuum_toroidal_field
    time_slice = equilibrium_ids.time_slice
    time_slice.resize(time_index + 1)
    time_slice = time_slice[0]
    global_quantities = time_slice.global_quantities
    boundary_outline = time_slice.boundary.outline
    profiles_1d = time_slice.profiles_1d
    profiles_2d = time_slice.profiles_2d
    profiles_2d.resize(profiles_2d_index + 1)
    profiles_2d = profiles_2d[profiles_2d_index]

    equilibrium_ids.time.resize(time_index + 1)
    equilibrium_ids.time = np.array(list(range(time_index + 1))) + 0.0
    equilibrium_ids.ids_properties.comment = "eqdsk python package"
    equilibrium_ids.ids_properties.name = eqdsk.name
    equilibrium_ids.ids_properties.homogeneous_time = 1
    global_quantities.ip = eqdsk.cplasma
    vacuum_toroidal_field.b0.resize(time_index + 1)
    vacuum_toroidal_field.b0[time_index] = np.array([eqdsk.bcentre] * (time_index + 1))
    vacuum_toroidal_field.r0 = eqdsk.xcentre
    global_quantities.psi_boundary = eqdsk.psibdry
    global_quantities.psi_axis = eqdsk.psimag
    global_quantities.psi_magnetic_axis = eqdsk.psimag
    global_quantities.magnetic_axis.r = eqdsk.xmag
    global_quantities.magnetic_axis.z = eqdsk.zmag
    profiles_1d.psi_norm = eqdsk.psinorm
    profiles_1d.psi = (eqdsk.psibdry - eqdsk.psimag) * eqdsk.psinorm + eqdsk.psimag
    profiles_1d.f_df_dpsi = eqdsk.ffprime
    profiles_1d.f = eqdsk.fpol
    profiles_1d.dpressure_dpsi = eqdsk.pprime
    profiles_1d.pressure = eqdsk.pressure
    profiles_1d.q = eqdsk.qpsi
    profiles_2d.psi = eqdsk.psi
    profiles_2d.grid.dim1 = np.linspace(
        eqdsk.xgrid1, eqdsk.xdim + eqdsk.xgrid1, eqdsk.nx
    )
    zmax = (2 * eqdsk.zmid + eqdsk.zdim) / 2
    profiles_2d.grid.dim2 = np.linspace(zmax - eqdsk.zdim, zmax, eqdsk.nz)
    boundary_outline.r = eqdsk.xbdry
    boundary_outline.z = eqdsk.zbdry

    db.put(equilibrium_ids)

    if eqdsk.nlim > 0:
        limiter_ids = ids_factory.wall()
        outline = limiter_ids.description_2d[0].limiter.unit[0].outline

        outline.r = eqdsk.xlim
        outline.z = eqdsk.zlim

        db.put(limiter_ids)

    if eqdsk.ncoil > 0:
        pf_active_ids = ids_factory.pf_active()

        for coil_id in range(eqdsk.ncoil):
            coil = pf_active_ids.coil[coil_id]

            coil.name = eqdsk.coil_names[coil_id]
            coil.current = eqdsk.Ic[coil_id]
            geometry = coil.geometry
            geometry.geometry_type = 2
            geometry.rectangle.width = eqdsk.xc[coil_id]
            geometry.rectangle.height = eqdsk.zc[coil_id]

        db.put(pf_active_ids)
