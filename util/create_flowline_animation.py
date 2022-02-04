#!/usr/bin/env python
# Copyright (C) 2021 Andy Aschwanden

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as colors

import numpy as np
import pylab as plt
import xarray as xr
import cf_xarray.units  # must be imported before pint_xarray
import pint_xarray
import os
import re
import ffmpeg

matplotlib.use("agg")

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.description = "Generating scripts for warming experiments."
parser.add_argument(
    "-b",
    "--basedir",
    dest="basedir",
    help="""Base directory.""",
    default="2021_12_fractures_melt",
)
parser.add_argument("FILE", nargs=1)
options = parser.parse_args()

infile = options.FILE[0]
basedir = options.basedir
smoothing_length = 1
m_id = re.search("id_(.+?)_", infile).group(1)

odir = os.path.join(basedir, "flowline", f"id_{m_id}")
if not os.path.isdir(odir):
    os.makedirs(odir)

ds = xr.open_dataset(infile)

time = ds["time"]
profile_axis_units = ds.variables["profile_axis"].attrs["units"]
ds = ds.pint.quantify({"profile_axis": profile_axis_units, "vonmises_stress": None})
ds = ds.pint.to(profile_axis="km").pint.dequantify()
config = ds.variables["pism_config"]
thickness_threshold = config.attrs["calving.thickness_calving.threshold"]

front_positions = {1985: 10.5, 2000: 12.5, 2003: 19.75, 2005: 24, 2010: 26}

for profile_id in ds.variables["profile_id"]:
    profile_name = ds.variables["profile_name"][profile_id]
    p_name = str(profile_name.values).lower().replace(" ", "_")
    profile_axis = ds.variables["profile_axis"][profile_id]
    x_min = profile_axis.min()
    x_max = profile_axis.max()

    print(f"Processing profile {profile_name.values}\n")

    thickness = ds.variables["thk"][profile_id, ...]
    surface = ds.variables["usurf"][profile_id, ...]
    surface = surface.where(thickness > thickness_threshold)
    topography = ds.variables["topg"][profile_id, ...]
    mask = ds.variables["mask"][profile_id, ...]
    dHdt = ds.variables["dHdt"][profile_id, ...]
    dHdt = dHdt.where(thickness > thickness_threshold)
    speed = ds.variables["velsurf_mag"][profile_id, ...]
    speed = speed.where(speed > 10)
    basal_melt = ds.variables["bmelt"][profile_id, ...]
    basal_melt = basal_melt.where((mask == 3) & (thickness > thickness_threshold))
    shelf_thickness = thickness.where((mask == 3) & (thickness > thickness_threshold))

    if smoothing_length > 1:
        basal_melt_rolling = (
            xr.DataArray(basal_melt).rolling(time=smoothing_length).mean()
        )
        dHdt_rolling = xr.DataArray(dHdt).rolling(time=smoothing_length).mean()
        speed_rolling = xr.DataArray(speed).rolling(time=smoothing_length).mean()
        shelf_thickness_rolling = (
            xr.DataArray(shelf_thickness).rolling(time=smoothing_length).mean()
        )
    else:
        basal_melt_rolling = basal_melt
        dHdt_rolling = dHdt
        speed_rolling = speed
        shelf_thickness_rolling = shelf_thickness
    bottom = surface - thickness

    cmap = plt.get_cmap("magma_r")
    cNorm = colors.Normalize(vmin=0, vmax=len(time))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)

    colors = cmap(np.arange(cmap.N))

    nt = len(time)
    for t, date in enumerate(time):
        print(f"""Processing {date.dt.strftime("%Y-%m-%d").values}""")
        fig, axs = plt.subplots(
            5,
            1,
            sharex="col",
            figsize=[10, 6],
            gridspec_kw=dict(height_ratios=[0.5, 2, 2, 2, 3]),
        )
        colorVals = scalarMap.to_rgba(range(0, t))
        axs[0].imshow([colors], extent=[profile_axis.min(), profile_axis.max(), 0, 1])
        pos = t / nt
        axs[0].axvline(
            pos * (profile_axis.max() - profile_axis.min()), color="k", lw=2.0
        )
        # axs[0].plot(profile_axis, speed_rolling[0:t, ...].T, color="k", lw=0.5)
        [
            axs[1].plot(
                profile_axis, speed_rolling[mt, ...].T, color=colorVals[mt], lw=0.5
            )
            for mt in range(0, t)
        ]
        axs[1].plot(profile_axis, speed_rolling[t, ...], color="k", lw=2.0)
        [
            axs[2].plot(
                profile_axis, basal_melt_rolling[mt, ...].T, color=colorVals[mt], lw=0.5
            )
            for mt in range(0, t)
        ]

        axs[2].plot(profile_axis, basal_melt_rolling[t, ...], color="k", lw=2.0)
        ax_dHdt = axs[2].twinx()
        ax_dHdt.plot(profile_axis, dHdt_rolling[t, ...], color="#238b45", lw=2.0)
        [
            axs[3].plot(
                profile_axis,
                shelf_thickness_rolling[mt, ...].T,
                color=colorVals[mt],
                lw=0.5,
            )
            for mt in range(0, t)
        ]

        axs[3].plot(profile_axis, shelf_thickness_rolling[t, ...], color="k", lw=2.0)

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
            color="k",
            linewidth=1.0,
        )
        axs[-1].fill_between(
            profile_axis,
            bottom[t, ...],
            surface[t, ...],
            color="#bdbdbd",
            linewidth=0.0,
        )

        axs[-1].fill_between(
            profile_axis,
            topography[t, ...] * 0 - 2000,
            topography[t, ...],
            color="#fdbe85",
            linewidth=0.3,
        )
        [
            axs[-1].axvline(dist, color=scalarMap.to_rgba(int((year - 1980) * 12)))
            for year, dist in front_positions.items()
        ]
        [
            axs[-1].text(
                dist, 1000, year, color=scalarMap.to_rgba(int((year - 1980) * 12))
            )
            for year, dist in front_positions.items()
        ]
        axs[-1].plot(profile_axis, topography[t, ...], color="k", lw=0.5)
        axs[-1].plot(profile_axis, bottom[t, ...], color="k", lw=0.5)
        axs[-1].plot(profile_axis, surface[t, ...], color="k", lw=0.5)
        axs[1].set_ylabel("""Speed\n(m/yr)""")
        axs[2].set_ylabel("""Subshelf\nmelt rate\n(m/yr)""")
        axs[3].set_ylabel("""Shelf thickness\n(m)""")
        axs[-1].set_xlabel(f"""Distance ({profile_axis.attrs["units"]})""")
        axs[-1].set_ylabel(f"""Altitude ({surface.attrs["units"]})""")
        ax_dHdt.set_ylabel("dHdt (m/yr)")
        axs[0].text(
            pos,
            1.2,
            date.dt.strftime("%Y-%m-%d").values,
            horizontalalignment="center",
            verticalalignment="center",
            transform=axs[0].transAxes,
        )
        axs[0].axes.xaxis.set_visible(False)
        axs[0].axes.yaxis.set_visible(False)
        axs[1].set_ylim(0, 15000)
        axs[2].set_ylim(100, 300)
        axs[3].set_ylim(0, 1000)
        axs[-1].set_ylim(-1500, 1000)
        ax_dHdt.set_ylim(-50, 50)
        for ax in axs:
            ax.set_xlim(x_min, x_max)
        fig.savefig(f"{odir}/{p_name}_{t:03d}.png", dpi=300)
        plt.close(plt.gcf())
        del fig, axs

    mdir = os.path.join(basedir, "flowline")
    (
        ffmpeg.input(f"{odir}/{p_name}_%03d.png")
        .output(
            f"{mdir}/{p_name}_id_{m_id}_hd1920.mp4",
            **{
                "s:v": "1920x1080",
                "c:v": "libx264",
                "crf": 20,
                "pix_fmt": "yuv420p",
                "framerate": 8,
                "r": 8,
            },
        )
        .run()
    )
ds.close()
