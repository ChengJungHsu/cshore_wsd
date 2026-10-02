import matplotlib.pyplot as plt
from plot_tools_cshore import *

parse_and_savenc()
model_A = 'OBPROF.nc'

x, zbi, zbf, zb_all, time, dates = extract_netcdf(model_A)
print(f'Simulation period  : {dates[0].strftime("%Y-%m-%d")}  to  {dates[-1].strftime("%Y-%m-%d")}')
print(f'Output timesteps   : {len(dates)}')
print(f'Cross-shore extent : {x[0]:.0f} to {x[-1]:.0f} m')

plot_profile_evolution(model_A)
plt.xlim(600,1000)
plt.savefig("evolution.png", dpi=300, bbox_inches="tight")
plot_bed_change(model_A)
plt.xlim(600,1000)
plt.savefig("dem_df.png", dpi=300, bbox_inches="tight")
plot_bed_change_point(model_A, 850)
plt.savefig("points_dem_df.png", dpi=300, bbox_inches="tight")
plot_profile_evolution_index(model_A,index=200)
plt.xlim(600,1000)
plot_profile_evolution_index(model_A,index=800)
plt.xlim(600,1000)
plot_profile_evolution_index(model_A,index=1400)
plt.xlim(600,1000)
plt.show()

# import datetime
# pos = dates.index(datetime.datetime(2019,1,1,0,0))
# with open("zb_20190101000000.txt", "w") as file:
    # # Loop through both variables simultaneously
    # for x, z in zip(x, zb_all[pos,:]):
        # # Write both variables separated by a tab and a newline character
        # file.write(f"{x}\t{z}\n")

