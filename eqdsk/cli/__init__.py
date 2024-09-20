# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Eqdsk CLI"""

import click
import numpy as np

from eqdsk.file import EQDSKInterface
from eqdsk.log import eqdsk_banner
from tests._helpers import read_strict_geqdsk  # noqa: PLC2701


@click.group()
@click.version_option()
def cli():
    """EQDSK cli tool.

    Contains useful utilities for display data in eqdsk files.
    """


@cli.command("show", no_args_is_help=True)
@click.argument("filepath", type=click.Path(exists=True))
@click.option(
    "--plot",
    default=False,
    is_flag=True,
    help="Display a plot the eqdsk plasma boundary.",
)
@click.option(
    "--plot-psi",
    default=False,
    is_flag=True,
    help="Display a plot the eqdsk psi data.",
)
def show(filepath, plot: bool, plot_psi: bool):  # noqa: FBT001
    """
    Reads an eqdsk and prints important parameters.

    May optionally plot the plasma boundary or psi data.

    matplotlib is required for plotting.
    """
    eq = EQDSKInterface.from_file(filepath, no_cocos=True)

    try:
        # tests if it is geqdsk readable/compatible
        read_strict_geqdsk(filepath)
        geqdsk_compatible = True
    except Exception:  # noqa: BLE001
        geqdsk_compatible = False

    print(eqdsk_banner())  # noqa: T201
    print(f"geqdsk compatible: {geqdsk_compatible}")  # noqa: T201
    print(eq)  # noqa: T201

    if plot or plot_psi:
        try:
            import matplotlib.pyplot as plt  # noqa: PLC0415
        except ImportError:
            print(  # noqa: T201
                "\n\nMatplotlib not installed, cannot plot.\n"
                "To resolve, try running:\n\n"
                "pip install matplotlib"
            )
            return

        fig, ax = plt.subplots()
        fig.set_size_inches(8, 10)

        ax.set_title(f"EQDSK: {eq.name}")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")

        # plot mag center
        ax.plot(eq.xmag, eq.zmag, color="red", marker="o")

        zoom = 1.2
        if plot_psi:
            # psi contour plot
            psi_x, psi_z = np.meshgrid(eq.x, eq.z)
            ax.contour(psi_x, psi_z, eq.psi, levels=20, colors="black")
            # TODO @oliverfunk: fix this legend, psi field not showen  # noqa: TD003
            ax.legend(["Magnetic center", "psi field"], loc="upper right")

            plot_max_x = eq.xmag + (eq.xdim / 2) * zoom
            plot_max_z = eq.zmag + (eq.zdim / 2) * zoom
        else:
            # plot lcfs
            ax.plot(eq.xbdry, eq.zbdry, color="blue", linestyle="-", linewidth=2)
            ax.legend(["Magnetic center", "Plasma Boundary"], loc="upper right")

            plot_max_x = max(eq.xbdry) * zoom
            plot_max_z = max(eq.zbdry) * zoom

        ax.set_xlim(0, plot_max_x)
        ax.set_ylim(-plot_max_z, plot_max_z)

        plt.show()
