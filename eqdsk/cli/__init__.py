# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Eqdsk CLI"""

import click
import numpy as np

from eqdsk.file import EQDSKInterface
from eqdsk.log import eqdsk_banner, eqdsk_warn
from tests._helpers import read_strict_geqdsk  # noqa: PLC2701


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

    ax.set_title(f"EQDSK:  {eq.name}")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")

    # plot mag center
    ax.plot(eq.xmag, eq.zmag, color="red", marker="o")

    return fig, ax, plt.show


@click.group()
@click.version_option()
def cli():
    """EQDSK cli

    A useful utility for displaying eqdsk file data.
    """


@cli.command("show", no_args_is_help=True)
@click.argument("filepath", type=click.Path(exists=True))
def show(filepath):
    """Reads and prints important parameters of the eqdsk."""
    eq = EQDSKInterface.from_file(filepath, no_cocos=True)

    try:
        # tests if it is geqdsk readable/compatible
        read_strict_geqdsk(filepath)
        geqdsk_compatible = True
    except Exception:  # noqa: BLE001
        geqdsk_compatible = False

    print(eqdsk_banner())  # noqa: T201
    print(f"geqdsk spec compliant: {geqdsk_compatible}")  # noqa: T201
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

    psi_x, psi_z = np.meshgrid(eq.x, eq.z)

    # TODO: should investigate the cause of this  # noqa: TD002, TD003
    if eq.psi.shape[0] != psi_x.shape[0]:
        eqdsk_warn(
            f"psi shape {eq.psi.shape} does match coords shape {psi_x.shape},"
            " attempting to swap"
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
