# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Eqdsk CLI"""

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Literal

import click
import numpy as np

from eqdsk.cocos import COCOS, KnownCOCOS
from eqdsk.file import IMAS_AVAIL, EQDSKInterface
from eqdsk.log import eqdsk_banner, eqdsk_warn

if IMAS_AVAIL:
    from eqdsk.file import DBEntry


def _imas_decoration(func: Callable[..., None]) -> Callable[..., None]:
    ti = click.option(
        "-ti", "--imas-time", "time", default=0, type=float, help="Time in"
    )
    tind = click.option(
        "-t-ind", "--imas-time-index", "t_ind", default=0, type=int, help="Time index"
    )
    pind = click.option(
        "-p-ind",
        "--imas-profile-2d-index",
        "p_ind",
        default=0,
        type=int,
        help="Profiles index",
    )

    dd = click.option(
        "-d", "--imas-dd-version", "dd_version", help="Data dictionary version"
    )

    return filepath(ti(tind(pind(dd(func)))))


class IMASPath(click.Path):
    """Path existence checking with imas avoidance"""

    def convert(self, value: str, *args, **kwargs):
        """Override Path check for imas databases"""  # noqa: DOC201
        if value.startswith("imas:"):
            return value
        return super().convert(value, *args, **kwargs)


filepath = click.argument("filepath_or_uri", type=IMASPath(exists=True))


def _setup_plotting(eq: EQDSKInterface):
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415
    except ImportError:
        raise ImportError(
            "Matplotlib not installed, cannot plot.\n"
            "To resolve, try running:\n"
            "pip install matplotlib"
        ) from None

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 8)

    fig.suptitle(f"EQDSK:  {eq.name}")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_aspect("equal")

    # plot mag center
    if None in {eq.xmag, eq.zmag}:
        eqdsk_warn("No xmag and zmag in eqdsk")
    else:
        ax.plot(eq.xmag, eq.zmag, color="red", marker="o", label="Magnetic centre")
        ax.legend(loc="upper right")

    return fig, ax, plt.show


@click.group()
@click.version_option()
def cli():
    """EQDSK cli

    Useful utilities for displaying eqdsk file data.
    """


def _imas_or_file_read(
    filepath_or_uri: str, tind: int, pind: int, time: float, dd_version: str
) -> EQDSKInterface:
    if filepath_or_uri.startswith("imas:") or filepath_or_uri.endswith(".nc"):
        with DBEntry(uri=filepath_or_uri, dd_version=dd_version, mode="r") as db:
            return EQDSKInterface.from_imas(
                db, time_index=tind, profiles_2d_index=pind, time=time, to_cocos=None
            )
    return EQDSKInterface.from_file(filepath_or_uri, no_cocos=True)


@cli.command("show", no_args_is_help=True)
@_imas_decoration
def show(filepath_or_uri: str, time: float, t_ind: int, p_ind: int, dd_version: str):
    """Reads and prints important parameters of the eqdsk.

    The below options only change how imas files are read.
    """
    eq = _imas_or_file_read(filepath_or_uri, t_ind, p_ind, time, dd_version)
    print(eqdsk_banner())  # noqa: T201
    print(eq)  # noqa: T201


@cli.command("plot", no_args_is_help=True)
@_imas_decoration
def plot_bdry(
    filepath_or_uri: str, t_ind: int, p_ind: int, time: float, dd_version: str
):
    """
    Plot the eqdsk plasma boundary.

    matplotlib is required for plotting.

    The below options only change how imas files are read.
    """
    eq = _imas_or_file_read(filepath_or_uri, t_ind, p_ind, time, dd_version)

    _fig, ax, show = _setup_plotting(eq)

    if 0 in {eq.xbdry.size, eq.zbdry.size}:
        eqdsk_warn("No xbdry and zbrdy in eqdsk")
    else:
        ax.plot(
            eq.xbdry,
            eq.zbdry,
            color="blue",
            linestyle="-",
            linewidth=2,
            label="Plasma Boundary",
        )

        zoom = 1.1
        plot_max_x = max(eq.xbdry) * zoom
        plot_max_z = max(eq.zbdry) * zoom

        ax.set_xlim(0, plot_max_x)
        ax.set_ylim(-plot_max_z, plot_max_z)

        ax.legend(loc="upper right")
    show()


@cli.command("plot-psi", no_args_is_help=True)
@_imas_decoration
def plot_psi(filepath_or_uri: str, t_ind: int, p_ind: int, time: float, dd_version: str):
    """
    Plot the eqdsk psi map.

    matplotlib is required for plotting.

    The below options only change how imas files are read.
    """
    eq = _imas_or_file_read(filepath_or_uri, t_ind, p_ind, time, dd_version)

    _fig, ax, show = _setup_plotting(eq)

    ax.set_title(
        "psi [V.s/rad (COCOS 1-8), V.s (COCOS 11-18)]",
        loc="right",
        fontdict={"fontsize": 10},
    )

    psi_x, psi_z = np.meshgrid(eq.x, eq.z, indexing="ij")

    if eq.psi.shape[0] != psi_x.shape[0]:
        eqdsk_warn(
            f"psi shape {eq.psi.shape} does match coords shape {psi_x.shape},"
            " attempting to swap..."
        )
        eq.psi = np.swapaxes(eq.psi, 0, 1)
    cont = ax.contour(psi_x, psi_z, eq.psi, levels=30, cmap="viridis")
    ax.clabel(cont, inline=True, fontsize=8, fmt="%d", colors="black")

    plot_max_x = None if eq.xmag is None else eq.xmag + (eq.xdim / 2)
    plot_max_z = None if eq.zmag is None else eq.zmag + (eq.zdim / 2)

    ax.set_xlim(0, plot_max_x)
    ax.set_ylim(None if plot_max_z is None else -plot_max_z, plot_max_z)

    show()


class COCOSOptionsHelp(click.Command):
    """COCOS options help formatter"""

    def format_help(self, ctx: click.Context, formatter: click.HelpFormatter):
        """Format the help with the Known COCOS options"""
        self.help = (self.help or "").format(
            "'" + "', '".join(KnownCOCOS.__members__.keys()) + "'"
        )
        super().format_help(ctx, formatter)

    def format_options(self, ctx: click.Context, formatter: click.HelpFormatter):
        """Writes all the options into the formatter if they exist."""
        opts, imas_opts = [], []

        for param in self.get_params(ctx):
            rv = param.get_help_record(ctx)
            if rv is not None:
                if "--imas" in rv[0]:
                    imas_opts.append(rv)
                else:
                    opts.append(rv)

        if opts:
            with formatter.section("Options"):
                formatter.write_dl(opts)

        if imas_opts:
            with formatter.section("IMAS Options"):
                formatter.write_dl(imas_opts)


def _dd_callback(ctx: click.Context, param, value: str | None) -> list[str] | None:  # noqa: ARG001
    return value.split(":") if isinstance(value, str) else value


@cli.command("convert", no_args_is_help=True, cls=COCOSOptionsHelp)
@filepath
@click.option(
    "-fmt",
    "--format",
    "format_",
    type=click.Choice(["json", "eqdsk", "imas"]),
    default="json",
    help="Format to save the eqdsk file in.",
)
@click.option("-f", "--from", "from_", help="COCOS format to read the eqdsk as.")
@click.option("-t", "--to", help="COCOS format to convert the eqdsk to.")
@click.option(
    "-q",
    "--qpsi-sign",
    type=click.Choice(["1", "-1"]),
    help="Sign of qpsi, "
    "required if the eqdsk has no qpsi data and you are converting between COCOS.",
)
@click.option("-u", "--imas-uri-out", "uri", help="File/uri output")
@click.option(
    "-m",
    "--imas-mode-out",
    "mode",
    type=click.Choice(["w", "a", "x"]),
    default="x",
    help="Output mode",
)
@click.option(
    "-d",
    "--imas-dd-version",
    "dd_version",
    callback=_dd_callback,
    help="I/O Data dictionary versions, I/O split with ':' eg '3.42.0:4.0.0'",
)
@click.option("-ti", "--imas-time", "time", default=0, type=float, help="Time in")
@click.option(
    "-t-ind",
    "--imas-time-index",
    "t_ind",
    default=(0, 0),
    type=(int, int),
    help="Time index in/out",
)
@click.option(
    "-p-ind",
    "--imas-profile-2d-index",
    "p_ind",
    default=(0, 0),
    type=(int, int),
    help="Profiles index in/out",
)
def convert(  # noqa: PLR0913, PLR0917
    filepath_or_uri: str,
    format_: Literal["eqdsk", "json", "imas"],
    from_: str | None,
    to: str | None,
    qpsi_sign: Literal["1", "-1"] | None,
    uri: str,
    mode: str,
    dd_version: tuple[str, ...] | None,
    time: float,
    t_ind: tuple[int, int],
    p_ind: tuple[int, int],
):
    """
    Conversion utilities for the eqdsk file.

    To save the file in a different format, use the -fmt (--format) option.

    To convert between COCOS versions, use the -f (--from) and -t (--to)
    options (both must be provided, unless the I/O is for imas).

    Valid values for --from and --to are:

      Integers: [1, 8] and [11, 18]

      Strings: {}

    The saved file will have '_out' suffixed to the filename.

    For non imas I/O:

        If -f and -t are not provided, the file will be read without COCOS.

        The specified "from" COCOS value must be valid for the eqdsk file.

        Use the `eqdsk show` command to see the valid COCOS's for file.

    For imas output --imas-uri-out can be specified to save to a specific database
    """  # noqa: DOC501
    if filepath_or_uri.startswith("imas:") or filepath_or_uri.endswith(".nc"):
        if from_:
            eqdsk_warn("from is not used as IMAS has a fixed COCOS")
        dv = dd_version[0] if isinstance(dd_version, Sequence) else None
        with DBEntry(uri=filepath_or_uri, mode="r", dd_version=dv) as db:
            eq = EQDSKInterface.from_imas(
                db, time_index=t_ind[0], profiles_2d_index=p_ind[0], time=time
            )
        if to is not None:
            if format_ == "imas":
                eqdsk_warn("Using IMAS db COCOS for output, ignoring '--to' argument")
            else:
                eq = eq.to_cocos(COCOS(to))
    elif from_ and to:
        # does validation of from and to values
        cc_fr = COCOS(from_)
        cc_to = COCOS(to)
        qsp = None if qpsi_sign is None else int(qpsi_sign) == 1
        eq = EQDSKInterface.from_file(
            filepath_or_uri,
            from_cocos=cc_fr.index,
            to_cocos=cc_to.index,
            qpsi_positive=qsp,
        )
    elif (from_ and to is None and format_ != "imas") or (from_ is None and to):
        # from_ imas check is not needed as it should be captured by first if
        raise click.BadParameter("Both --from and --to must be provided")

    else:
        eq = EQDSKInterface.from_file(filepath_or_uri, no_cocos=True)

    fp = Path(filepath_or_uri.replace("imas:hdf5?path=", ""))
    output_path = fp.with_stem(f"{fp.stem}_out")
    if format_ == "imas":
        if from_ is None and not (
            filepath_or_uri.startswith("imas:") or filepath_or_uri.endswith(".nc")
        ):
            raise click.BadParameter("--from must be provided for conversion to imas")
        if eq._cocos is None:  # noqa: SLF001
            eq.identify(as_cocos=from_)

        if uri is None:
            uri = output_path.with_suffix(".nc").as_posix()
        dv = dd_version[-1] if isinstance(dd_version, Sequence) else None
        with DBEntry(uri=uri, mode=mode, dd_version=dv) as db:
            eq.write(
                db,
                file_format=format_,
                imas_kwargs={
                    "time_index": t_ind[1],
                    "profiles_2d_index": p_ind[1],
                },
            )
    else:
        eq.write(output_path.as_posix(), file_format=format_)
