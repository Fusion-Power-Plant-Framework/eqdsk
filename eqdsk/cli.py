# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Eqdsk CLI"""

from pathlib import Path

import click
import numpy as np

from eqdsk.cocos import COCOS, KnownCOCOS
from eqdsk.file import EQDSKInterface
from eqdsk.log import eqdsk_banner, eqdsk_warn


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
    ax.plot(eq.xmag, eq.zmag, color="red", marker="o")

    return fig, ax, plt.show


@click.group()
@click.version_option()
def cli():
    """EQDSK cli

    Useful utilities for displaying eqdsk file data.
    """


@cli.command("show", no_args_is_help=True)
@click.argument("filepath", type=click.Path(exists=True))
def show(filepath):
    """Reads and prints important parameters of the eqdsk."""
    eq = EQDSKInterface.from_file(filepath, no_cocos=True)

    print(eqdsk_banner())  # noqa: T201
    print(eq)  # noqa: T201


@cli.command("plot", no_args_is_help=True)
@click.argument("filepath", type=click.Path(exists=True))
def plot_bdry(filepath):
    """
    Plot the eqdsk plasma boundary.

    matplotlib is required for plotting.
    """
    eq = EQDSKInterface.from_file(filepath, no_cocos=True)

    _fig, ax, show = _setup_plotting(eq)

    ax.plot(eq.xbdry, eq.zbdry, color="blue", linestyle="-", linewidth=2)
    ax.legend(["Magnetic center", "Plasma Boundary"], loc="upper right")

    zoom = 1.1
    plot_max_x = max(eq.xbdry) * zoom
    plot_max_z = max(eq.zbdry) * zoom

    ax.set_xlim(0, plot_max_x)
    ax.set_ylim(-plot_max_z, plot_max_z)

    show()


@cli.command("plot-psi", no_args_is_help=True)
@click.argument("filepath", type=click.Path(exists=True))
def plot_psi(filepath):
    """
    Plot the eqdsk psi map.

    matplotlib is required for plotting.
    """
    eq = EQDSKInterface.from_file(filepath, no_cocos=True)

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

    ax.legend(["Magnetic center"], loc="upper right")

    plot_max_x = eq.xmag + (eq.xdim / 2)
    plot_max_z = eq.zmag + (eq.zdim / 2)

    ax.set_xlim(0, plot_max_x)
    ax.set_ylim(-plot_max_z, plot_max_z)

    show()


class COCOSOptionsHelp(click.Command):
    """COCOS options help formatter"""

    def format_help(self, ctx, formatter):
        """Format the help with the Known COCOS options"""
        self.help = self.help.format(
            "'" + "', '".join(KnownCOCOS.__members__.keys()) + "'"
        )
        super().format_help(ctx, formatter)


@cli.command("convert", no_args_is_help=True, cls=COCOSOptionsHelp)
@click.argument("filepath", type=click.Path(exists=True))
@click.option(
    "-fmt",
    "--format",
    type=click.Choice(["json", "eqdsk"]),
    default="json",
    help="Format to save the eqdsk file in.",
)
@click.option(
    "-f",
    "--from",
    "from_",
    help="COCOS format to read the eqdsk as.",
)
@click.option(
    "-t",
    "--to",
    help="COCOS format to convert the eqdsk to.",
)
@click.option(
    "-q",
    "--qpsi-sign",
    type=click.Choice(["1", "-1"]),
    help="Sign of qpsi, "
    "required if the eqdsk has no qpsi data and you are converting between COCOS.",
)
def convert(
    filepath: str,
    format: str,  # noqa: A002
    from_: str | None,
    to: str | None,
    qpsi_sign: str | None,
):
    """
    Conversion utilities for the eqdsk file.

    To save the file in a different format, use the -fmt (--format) option.

    To convert between COCOS versions, use the -f (--from) and -t (--to)
    options (both must be provided).

    If -f and -t are not provided, the file will be read without COCOS.

    Valid values for --from and --to are:

      Integers: [1, 8] and [11, 18]

      Strings: {}

    The specified "from" COCOS value must be valid for the eqdsk file.

    Use the `eqdsk show` command to see the valid COCOS's for file.

    The saved file will have _out suffixed to the filename.
    """  # noqa: DOC501
    if from_ and to:
        # does validation of from and to values
        cc_fr = COCOS(from_)
        cc_to = COCOS(to)
        qsp = None if qpsi_sign is None else int(qpsi_sign) == 1
        eq = EQDSKInterface.from_file(
            filepath,
            from_cocos=cc_fr.index,
            to_cocos=cc_to.index,
            qpsi_positive=qsp,
        )
    elif (from_ and to is None) or (from_ is None and to):
        raise click.BadParameter("Both --from and --to must be provided")
    else:
        eq = EQDSKInterface.from_file(filepath, no_cocos=True)

    output_path = Path(filepath).with_stem(f"{Path(filepath).stem}_out").as_posix()
    eq.write(output_path, file_format=format)
