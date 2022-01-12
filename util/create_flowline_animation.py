#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
import pylab as plt
import xarray as xr
import cf_xarray.units  # must be imported before pint_xarray
import pint_xarray
from pint_xarray import unit_registry as ureg
import os
import re
import ffmpeg

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Generating scripts for warming experiments."
parser.add_argument(
    "-b",
    "--basedir",
    dest="basedir",
    help="""Base directory.""",
    default="2021_12_all",
)
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()

infile = options.FILE[0]
basedir = options.basedir
smoothing_length = 13

m_id = re.search("id_(.+?)_", infile).group(1)

odir = os.path.join(basedir, "flowline", f"id_{m_id}")
if not os.path.isdir(odir):
    os.makedirs(odir)

ds = xr.open_dataset(infile)

time = ds["time"].dt.strftime("%Y-%m-%d").values
profile_axis_units = ds.variables["profile_axis"].attrs["units"]
ds = ds.pint.quantify({"profile_axis": profile_axis_units, "vonmises_stress": None})
ds = ds.pint.to(profile_axis="km").pint.dequantify()


for profile_id in ds.variables["profile_id"]:
    profile_name = ds.variables["profile_name"][profile_id]
    profile_axis = ds.variables["profile_axis"][profile_id]
    x_min = profile_axis.min()
    x_max = profile_axis.max()
    print(f"Processing profile {profile_name.values}\n")

    surface = ds.variables["usurf"][profile_id, ...]
    topography = ds.variables["topg"][profile_id, ...]
    thickness = ds.variables["thk"][profile_id, ...]
    dHdt = ds.variables["dHdt"][profile_id, ...]
    dHdt = dHdt.where(thickness > 10)
    speed = ds.variables["velsurf_mag"][profile_id, ...]
    speed = speed.where(speed > 10)
    basal_melt = ds.variables["bmelt"][profile_id, ...]
    basal_melt = basal_melt.where(basal_melt > 10)

    speed_rolling = xr.DataArray(speed).rolling(time=smoothing_length).mean()
    basal_melt_rolling = xr.DataArray(basal_melt).rolling(time=smoothing_length).mean()
    dHdt_rolling = xr.DataArray(dHdt).rolling(time=smoothing_length).mean()

    bottom = surface - thickness
    for t, date in enumerate(time):
        print(f"Processing {date}")
        fig, axs = plt.subplots(
            3,
            1,
            sharex="col",
            figsize=[6.2, 6],
            gridspec_kw=dict(height_ratios=[2, 2, 3]),
        )
        axs[0].plot(profile_axis, speed_rolling[0:t, ...].T, color="0.75", lw=0.5)
        axs[0].plot(profile_axis, speed_rolling[t, ...], color="k", lw=2.0)
        axs[1].plot(profile_axis, basal_melt_rolling[0:t, ...].T, color="0.75", lw=0.5)
        axs[1].plot(profile_axis, basal_melt_rolling[t, ...], color="k", lw=2.0)
        ax_dHdt = axs[1].twinx()
        ax_dHdt.plot(profile_axis, dHdt_rolling[t, ...], color="b", lw=2.0)

        axs[-1].fill_between(
            profile_axis,
            topography[t, ...] * 0 - 2000,
            topography[t, ...] * 0,
            color="#c6dbef",
            linewidth=0.3,
        )

        axs[-1].fill_between(
            profile_axis,
            bottom[t, ...],
            surface[t, ...],
            color="#bdbdbd",
            linewidth=0.5,
        )
        axs[-1].fill_between(
            profile_axis,
            topography[t, ...] * 0 - 2000,
            topography[t, ...],
            color="#fdbe85",
            linewidth=0.3,
        )
        axs[-1].plot(profile_axis, topography[t, ...], color="k", lw=0.5)
        axs[-1].plot(profile_axis, bottom[t, ...], color="k", lw=0.5)
        axs[-1].plot(profile_axis, surface[t, ...], color="k", lw=0.5)
        axs[0].set_ylabel(f"""Speed\n({speed.attrs["units"]})""")
        axs[1].set_ylabel(f"""Subshelf melt rate\n({basal_melt.attrs["units"]})""")

        axs[-1].set_xlabel(f"""Distance ({profile_axis.attrs["units"]})""")
        axs[-1].set_ylabel(f"""Altitude ({surface.attrs["units"]})""")

        axs[0].text(
            0.0,
            1.1,
            date,
            horizontalalignment="left",
            verticalalignment="center",
            transform=axs[0].transAxes,
        )
        axs[0].set_ylim(0, 20000)
        axs[1].set_ylim(100, 300)
        axs[-1].set_ylim(-1500, 1000)
        ax_dHdt.set_ylim(-50, 50)
        for ax in axs:
            ax.set_xlim(x_min, x_max)
        fig.savefig(f"{odir}/flowline_{t:03d}.png", dpi=300)
        plt.close(plt.gcf())
        del fig
    ds.close()

    (
        ffmpeg.input(f"{odir}/flowline_%03d.png")
        .output(
            f"{odir}/jib_flowline_id_{m_id}_hd1920.mp4",
            **{
                "s:v": "1920x1080",
                "c:v": "libx264",
                "crf": 20,
                "pix_fmt": "yuv420p",
                "framerate": 10,
                "r": 10,
            },
        )
        .run()
    )
