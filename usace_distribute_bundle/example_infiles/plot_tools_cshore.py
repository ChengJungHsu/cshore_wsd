from datetime import datetime, timedelta
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


REF_DATE = datetime(2019, 1, 1, 0, 0)
MHW = 1.08  # Mean High Water elevation (m)

# Standard colors for Models A, B, C, D
MODEL_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
def parse_and_savenc(file_path="OBPROF"):
    time_steps = []
    times = []
    zb_all = []
    with open(file_path, "r") as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        # Parse header line: contains "flag/ID", "number of rows", "time in seconds"
        parts = line.split()
        if len(parts) == 3:
            try:
                num_rows = int(parts[1])
                time_sec = float(parts[2])
            except ValueError:
                i += 1
                continue

            x_vals = []
            z_vals = []

            # Read the specified number of data rows following the header
            for j in range(1, num_rows + 1):
                if i + j < len(lines):
                    data_parts = lines[i + j].split()
                    if len(data_parts) >= 2:
                        x_vals.append(float(data_parts[0]))
                        z_vals.append(float(data_parts[1]))

            time_steps.append({"time": time_sec, "x": x_vals, "z": z_vals})
            times.append(time_sec)
            zb_all.append(z_vals)
            # Move index past the read rows
            i += num_rows + 1
        else:
            i += 1

    x = np.array(time_steps[0]["x"])
    times = np.array(times)
    zb_all = np.array(zb_all)
    
    times = times[::6]
    zb_all =zb_all[::6]

    # Save to NetCDF file at file_path
    with nc.Dataset(file_path+".nc", "w", format="NETCDF4") as ds:
        # Create dimensions: 'time' (len(times)) and 'x' (len(x))
        ds.createDimension("time", len(times))
        ds.createDimension("x", len(x))

        # Create coordinate variables
        var_times = ds.createVariable("times", "f8", ("time",))
        var_x = ds.createVariable("x", "f8", ("x",))

        # Create main 2D variable (shape: len(times) x len(x))
        var_zb = ds.createVariable("zb_all", "f8", ("time", "x"))

        # Add optional variable attributes
        var_times.units = "seconds"
        var_x.units = "meters"
        var_zb.units = "meters"

        # Write data array values
        var_times[:] = times
        var_x[:] = x
        var_zb[:, :] = zb_all

    return times, x, zb_all
    
def extract_netcdf(fname):
    """Load cross-shore profile data from an AeoLiS netCDF output file.

    Parameters
    ----------
    fname : str
        Path to netCDF output file

    Returns
    -------
    x : ndarray (nx,)
        Cross-shore coordinates
    zbi : ndarray (nx,)
        Initial bed level
    zbf : ndarray (nx,)
        Final bed level
    zb_all : ndarray (nt, nx)
        Bed level at all output times
    rhoveg_all : ndarray (nt, nx) or None
        Vegetation density at all output times (None if not in file)
    time : ndarray (nt,)
        Time in seconds from reference date
    dates : list of datetime
        Datetime objects for each output step
    """
    dataset = nc.Dataset(fname, "r")
    # Extract variables as NumPy arrays
    times = dataset.variables["times"][:]
    x = dataset.variables["x"][:]
    zb_all = dataset.variables["zb_all"][:]
    
    zbi = zb_all[0, :]
    zbf = zb_all[-1, :]
    dates = [REF_DATE + timedelta(seconds=float(t)) for t in times]
    return x, zbi, zbf, zb_all, times, dates


def _add_colorbar(fig, ax, dates):
    """Add a viridis time colorbar with year tick labels."""
    norm = mpl.colors.Normalize(vmin=0, vmax=len(dates) - 1)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=plt.cm.viridis)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label('Dates', fontsize=12, fontweight='bold')
    tick_locs = np.linspace(0, len(dates) - 1, min(6, len(dates))).astype(int)
    cbar.set_ticks(tick_locs)
    cbar.set_ticklabels([str(dates[i].year) for i in tick_locs])
    return cbar


def plot_profile_evolution(fname, mhw=MHW):
    """Plot cross-shore bed profile evolution colored by time (viridis).

    Each profile line is colored by the viridis colormap progressing from
    early (dark purple) to late (yellow).

    Parameters
    ----------
    fname : str
        Path to netCDF output file
    mhw : float
        Mean high water elevation (m) shown as a grey dashed reference line
    """
    x, _, _, zb_all, time, dates = extract_netcdf(fname)
    colors_arr = plt.cm.viridis(np.linspace(0, 1, len(time)))

    fig, ax = plt.subplots(figsize=(8, 5))
    for i in range(len(time)):
        ax.plot(x, zb_all[i, :], color=colors_arr[i])
    ax.axhline(y=mhw, color='grey', linestyle='--', label='Mean High Water')
    ax.set_ylabel('Elevation (m)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Cross-shore Distance (m)', fontsize=13, fontweight='bold')
    ax.grid()
    ax.legend(loc='upper left')
    _add_colorbar(fig, ax, dates)
    plt.tight_layout()
    return fig, ax
    
def plot_profile_evolution_index(fname, mhw=MHW,index=5000):
    """Plot cross-shore bed profile evolution colored by time (viridis).

    Each profile line is colored by the viridis colormap progressing from
    early (dark purple) to late (yellow).

    Parameters
    ----------
    fname : str
        Path to netCDF output file
    mhw : float
        Mean high water elevation (m) shown as a grey dashed reference line
    """
    x, _, _, zb_all, time, dates = extract_netcdf(fname)
    zb_all=zb_all[0:index,:]
    time=time[0:index]
    dates=dates[0:index]
    colors_arr = plt.cm.viridis(np.linspace(0, 1, len(time)))

    fig, ax = plt.subplots(figsize=(8, 5))
    for i in range(len(time)):
        ax.plot(x, zb_all[i, :], color=colors_arr[i])
    ax.axhline(y=mhw, color='grey', linestyle='--', label='Mean High Water')
    ax.set_ylabel('Elevation (m)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Cross-shore Distance (m)', fontsize=13, fontweight='bold')
    ax.set_title(f"{dates[-1]}")
    ax.grid()
    ax.legend(loc='upper left')
    _add_colorbar(fig, ax, dates)
    plt.tight_layout()
    return fig, ax    

def plot_bed_change_point(fname,index=500):
    """Plot net bed-level change (final minus initial) with shading.

    Blue shading = deposition, red shading = erosion.

    Parameters
    ----------
    fname : str
        Path to netCDF output file
    """
    x, zbi, zbf, zb_all, times, dates = extract_netcdf(fname)
    dz = zb_all[:,index] - zbi[index]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(dates, dz, label=f"x = {x[index]}m" ,color='steelblue', linewidth=2)
    ax.set_ylabel('Bed-level Change (m)', fontsize=13, fontweight='bold')
    ax.set_xlabel('time', fontsize=13, fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid()
    plt.tight_layout()
    return fig, ax

def plot_bed_change(fname):
    """Plot net bed-level change (final minus initial) with shading.

    Blue shading = deposition, red shading = erosion.

    Parameters
    ----------
    fname : str
        Path to netCDF output file
    """
    x, zbi, zbf, _, _, _ = extract_netcdf(fname)
    dz = zbf - zbi

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, dz, color='steelblue', linewidth=2)
    ax.axhline(0, color='k', linewidth=0.8, linestyle='--')
    ax.fill_between(x, dz, 0, where=(dz > 0),
                    color='steelblue', alpha=0.3, label='Deposition')
    ax.fill_between(x, dz, 0, where=(dz < 0),
                    color='salmon', alpha=0.5, label='Erosion')
    ax.set_ylabel('Bed-level Change (m)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Cross-shore Distance (m)', fontsize=13, fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid()
    plt.tight_layout()
    return fig, ax


def _circular_mean_deg(angles_deg, axis=0):
    """Compute the circular mean of angles in degrees."""
    radians = np.deg2rad(angles_deg)
    sin_mean = np.nanmean(np.sin(radians), axis=axis)
    cos_mean = np.nanmean(np.cos(radians), axis=axis)
    mean_rad = np.arctan2(sin_mean, cos_mean)
    return np.mod(np.rad2deg(mean_rad), 360)


def plot_windrose_and_timeseries(fname):
    """Plot a windrose and wind speed/direction timeseries from AeoLiS output.

    Parameters
    ----------
    fname : str
        Path to netCDF output file
    """
    _, _, _, _, _, dates = extract_netcdf(fname)

    ds = nc.Dataset(fname, 'r')
    if 'uw' not in ds.variables:
        print(f'No wind speed variable "uw" found in {os.path.basename(fname)}')
        ds.close()
        return None, None

    if 'udir' not in ds.variables:
        print(f'No wind direction variable "udir" found in {os.path.basename(fname)}')
        ds.close()
        return None, None

    uw = ds['uw'][:]
    ud = ds['udir'][:]
    ds.close()

    if uw.ndim > 1:
        uw_flat = uw.reshape(uw.shape[0], -1)
    else:
        uw_flat = uw
    if ud.ndim > 1:
        ud_flat = ud.reshape(ud.shape[0], -1)
    else:
        ud_flat = ud

    mean_speed = np.nanmean(uw_flat, axis=1)
    mean_direction = _circular_mean_deg(ud_flat, axis=1)

    theta = np.deg2rad(ud_flat.flatten())
    speed_weights = uw_flat.flatten()
    valid = ~np.isnan(theta) & ~np.isnan(speed_weights)
    theta = theta[valid]
    speed_weights = speed_weights[valid]

    if theta.size == 0:
        print(f'No valid wind data found in {os.path.basename(fname)}')
        return None, None

    bins = np.linspace(0, 2 * np.pi, 17)
    counts, _ = np.histogram(theta, bins=bins, weights=speed_weights)
    widths = np.diff(bins)
    centers = bins[:-1] + 0.5 * widths

    fig = plt.figure(figsize=(14, 5))
    ax0 = fig.add_subplot(1, 2, 1, projection='polar')
    ax0.bar(centers, counts, width=widths, bottom=0.0, align='center', edgecolor='k', alpha=0.7)
    ax0.set_theta_zero_location('N')
    ax0.set_theta_direction(-1)
    ax0.set_xticks(np.deg2rad(np.arange(0, 360, 45)))
    ax0.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
    ax0.set_title('Windrose', fontsize=13, fontweight='bold')

    ax1 = fig.add_subplot(1, 2, 2)
    ax1.plot(dates, mean_speed, color='steelblue', linewidth=2, label='Mean Speed (m/s)')
    ax1.set_xlabel('Date', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean Speed (m/s)', color='steelblue', fontsize=12, fontweight='bold')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(dates, mean_direction, color='darkorange', linewidth=2, label='Mean Direction (deg)')
    ax2.set_ylabel('Mean Direction (deg)', color='darkorange', fontsize=12, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='darkorange')

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc='upper left')

    fig.autofmt_xdate()
    plt.tight_layout()
    return fig, (ax0, ax1)


def plot_dune_growth(fname, dune_toe_elevation=4.4):
    """Plot the cross-shore location of the dune toe elevation over time.

    Parameters
    ----------
    fname : str
        Path to netCDF output file
    dune_toe_elevation : float, optional
        Elevation threshold for the dune toe. If None, uses the mean high
        water elevation constant `MHW`.
    """
    x, _, _, zb_all, _, dates = extract_netcdf(fname)

    toe_x = np.full(len(dates), np.nan)
    for i in range(len(dates)):
        profile = zb_all[i, :]
        ix = np.where(profile >= dune_toe_elevation)[0]
        if ix.size > 0:
            j = ix[0]
            if j == 0:
                toe_x[i] = x[0]
            else:
                x0, x1 = x[j - 1], x[j]
                z0, z1 = profile[j - 1], profile[j]
                if z1 != z0:
                    try:
                        toe_x[i] = x0 + (dune_toe_elevation - z0) * (x1 - x0) / (z1 - z0)
                    except:
                        toe_x[i] = np.nan
                else:
                    toe_x[i] = x1

    years = [(d - dates[0]).days / 365.25 for d in dates]
    fig, ax = plt.subplots(figsize=(5, 6))
    ax.plot(toe_x, years, color='steelblue', linewidth=2)
    ax.set_ylabel('Time (years)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Dune Toe \nCross-shore Location (m)', fontsize=13, fontweight='bold')
    ax.grid()
    plt.tight_layout()
    return fig, ax

# def parse_and_plot_bathymetry(file_path):
    # time_steps = []
    # times = []
    # zb_all = []
    # with open(file_path, 'r') as file:
        # lines = file.readlines()
        
    # i = 0
    # while i < len(lines):
        # line = lines[i].strip()
        # if not line:
            # i += 1
            # continue
        
        # # Parse header line: contains "flag/ID", "number of rows", "time in seconds"
        # parts = line.split()
        # if len(parts) == 3:
            # try:
                # num_rows = int(parts[1])
                # time_sec = float(parts[2])
            # except ValueError:
                # i += 1
                # continue
            
            # x_vals = []
            # z_vals = []
            
            # # Read the specified number of data rows following the header
            # for j in range(1, num_rows + 1):
                # if i + j < len(lines):
                    # data_parts = lines[i + j].split()
                    # if len(data_parts) >= 2:
                        # x_vals.append(float(data_parts[0]))
                        # z_vals.append(float(data_parts[1]))
            
            # time_steps.append({
                # 'time': time_sec,
                # 'x': x_vals,
                # 'z': z_vals
            # })
            # times.append(time_sec)
            # zb_all.append(z_vals)
            # # Move index past the read rows
            # i += num_rows + 1
        # else:
            # i += 1
    # x = np.array(time_steps[0]['x'])
    # times = np.array(times)
    # zb_all = np.array(zb_all)
    # return times, x, zb_all


# def plot_vegetation_evolution(fname):
    # """Plot vegetation density evolution colored by time (viridis).

    # Parameters
    # ----------
    # fname : str
        # Path to netCDF output file
    # """
    #>>> x, _, _, _, rhoveg_all, time, dates = extract_netcdf(fname)
    # if rhoveg_all is None:
        # print(f'No vegetation output found in {os.path.basename(fname)}')
        # return None, None
    # veg_x = np.full(len(dates), np.nan)
    # for i in range(len(dates)):
        # veg = rhoveg_all[i, :]
        # target = 0.2
        # ix = np.where(veg >= target)[0]
        # if ix.size > 0:
            # j = ix[0]
            # if j == 0:
                # veg_x[i] = x[0]
            # else:
                # x0, x1 = x[j - 1], x[j]
                # z0, z1 = veg[j - 1], veg[j]
                # if z1 != z0:
                    # veg_x[i] = x0 + (target - z0) * (x1 - x0) / (z1 - z0)
                # else:
                    # veg_x[i] = x1

    # # colors_arr = plt.cm.viridis(np.linspace(0, 1, len(time)))
    # years = [(d - dates[0]).days / 365.25 for d in dates]

    # fig, ax = plt.subplots(figsize=(5, 6))
    # # for i in range(len(time)):
    # ax.plot(veg_x, years) #, color=colors_arr[i])
    # ax.set_ylabel('Time', fontsize=13, fontweight='bold')
    # ax.set_xlabel('First Cross-Shore \n Location of Established Vegetation (m)', fontsize=13, fontweight='bold')
    # # ax.set_ylim([0, 1])
    # ax.grid()
    # # _add_colorbar(fig, ax, dates)
    # plt.tight_layout()
    # return fig, ax
