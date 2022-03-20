import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import netCDF4 as nc

import PyFVCOM as pf

import flux_calc_funcs as fcf

bathy_limit = 10
static_depth = False
salinity_below_lim = -4
salinity_above_lim = 6
time_in_lim = 0.8


#rect_area = [[-4.25, 50.12] , [-4.05,50.5]]
rect_area = [[-4.25, 50], [-4, 50.25]]
area_ident = 'landsch'

#mod_filestr = '/home/mike/Experiments/Rich_paper/Data/run2_mod/tamar_v2_0001.nc'
mod_filestr = '/data/sthenno1/backup/mbe/Models/FVCOM/tamar_v2/output/phys_only_2017/2016/07/tamar_v2_0001.nc'

cnst_wind = 7.1

##################################################################
# Get the buoy wind, temp, sal
this_year = 2016
#buoy_file = '/home/mike/Experiments/Rich_paper/Data/l4_cont_data_{}.txt'.format(this_year)
buoy_file = '/data/sthenno1/backup/mbe/Data/WCO_data/L4/Buoy_data/l4_cont_data_{}.txt'.format(this_year)

buoy_raw = np.loadtxt(buoy_file, skiprows=1, dtype=str)

year = np.asarray(buoy_raw[:,0], dtype=int)
jd = np.asarray(buoy_raw[:,1], dtype=int)
time_raw = np.asarray(buoy_raw[:,2], dtype=float)

hour = np.floor(time_raw)
minute = np.floor(100*(time_raw-hour))

buoy_dt =[]
for i in np.arange(0, len(year)):
    this_str = '{}{}_{:02d}:{:02d}'.format(year[i],jd[i],int(hour[i]),int(minute[i]))
    buoy_dt.append(dt.datetime.strptime(this_str, '%Y%j_%H:%M'))
buoy_dt = np.asarray(buoy_dt)

wind_u = np.asarray(buoy_raw[:,-8], dtype='float')
wind_v = np.asarray(buoy_raw[:,-6], dtype='float')

remove_duff = wind_u == -999.99

#buoy_dt = buoy_dt[~remove_duff]
#wind_spd = np.sqrt(wind_u[~remove_duff]**2 + wind_v[~remove_duff]**2)

wrf_l4_data = np.load('wrf_l4_wnd.npy', allow_pickle=True)
buoy_dt = wrf_l4_data[0]
buoy_dt, buoy_ind = np.unique(buoy_dt, return_index=True)

wind_spd = np.asarray(wrf_l4_data[1], dtype=float)[buoy_ind]

"""
buoy_temp = np.asarray(buoy_raw[:,11], dtype='float')[~remove_duff]
buoy_sal = np.asarray(buoy_raw[:,13], dtype='float')[~remove_duff]

remove_ts = np.logical_and(buoy_temp < 0, buoy_sal < 25)
buoy_temp = buoy_temp[~remove_ts]
buoy_sal = buoy_sal[~remove_ts]
buoy_dt_ts = buoy_dt[~remove_ts]
"""
del buoy_raw


###############################################################
# construct interpolated PC02 etc

#data_dir = '/home/mike/Experiments/Rich_paper/Data/'
data_dir = '/data/sthenno1/backup/mbe/Experiments/Rich_paper/Data/'
l4_fC02_file = '{}Quest_2016_fCO2_L4_only.csv'.format(data_dir)

# Get the C02 measurements

fC02_data = np.loadtxt(l4_fC02_file, delimiter=',', skiprows=1)
fC02_dt_bst = np.asarray([dt.datetime(int(this_row[0]), int(this_row[1]), int(this_row[2]), int(this_row[3]), int(this_row[4]), int(this_row[5])) for this_row in fC02_data])
fC02_dt_raw = np.asarray([this_dt - dt.timedelta(hours=1) for this_dt in fC02_dt_bst])
fC02_day = np.asarray([dt.datetime(2016, this_dt.month, this_dt.day) for this_dt in fC02_dt_raw])
fC02_day_u = np.unique(fC02_day)

fC02_dt = []
fC02 = []
pC02 = []
pC02_atmos = []
fC02_lon = []
fC02_lat = []

for this_day in fC02_day_u:
    choose_dt = fC02_day == this_day
    
    poss_dates = fC02_dt_raw[choose_dt]
    fC02_dt.append(poss_dates[int(np.floor(len(poss_dates)/2))])

    fC02.append(np.mean(fC02_data[choose_dt,8]))
    pC02.append(np.mean(fC02_data[choose_dt,6]))
    pC02_atmos.append(np.mean(fC02_data[choose_dt,7]))

    fC02_lon.append(np.mean(fC02_data[choose_dt,13]))
    fC02_lat.append(np.mean(fC02_data[choose_dt,12]))

fC02_dt = np.asarray(fC02_dt)
fC02 = np.asarray(fC02)
pC02 = np.asarray(pC02)
pC02_atmos = np.asarray(pC02_atmos)
fC02_lon = np.asarray(fC02_lon)
fC02_lat = np.asarray(fC02_lat)

###############################################################
# Get the model temp, sal, zeta series
#tamar_data_dir = '/home/mike/Experiments/Rich_paper/Data/temp_sal_only'
tamar_data_dir = '/data/sthenno1/backup/mbe/Models/FVCOM/tamar_v2/output/phys_only_2017/2016'
tamar_data = nc.Dataset('{}/temp.nc'.format(tamar_data_dir))

tamar_grid = pf.read.FileReader(mod_filestr)

box_valid = np.logical_and(np.logical_and(tamar_grid.grid.lon >= rect_area[0][0],
		tamar_grid.grid.lon <= rect_area[1][0]), np.logical_and(tamar_grid.grid.lat >= rect_area[0][1],
		tamar_grid.grid.lat <= rect_area[1][1]))

mod_datestr = np.asarray([b''.join(this_str).decode('utf-8') for this_str in tamar_data.variables['Times'][:]])
mod_dt = np.asarray([dt.datetime.strptime(this_str, '%Y-%m-%dT%H:%M:%S.000000') for this_str in mod_datestr])
mod_sal = np.squeeze(tamar_data['salinity'][:])
mod_temp = np.squeeze(tamar_data['temp'][:])
mod_zeta = tamar_data['zeta'][:]

tamar_data.close()

# There are a few double entries from the concatenation, remove these
remove_rows = np.where(mod_dt[0:-1] == mod_dt[1:])[0]

mod_dt = np.delete(mod_dt, remove_rows)
mod_sal = np.delete(mod_sal, remove_rows, axis=0)
mod_temp = np.delete(mod_temp, remove_rows, axis=0)
mod_zeta = np.delete(mod_zeta, remove_rows, axis=0)

# Save the L4 sal and temp for use elsewhere
l4_mean = np.squeeze(tamar_grid.closest_node([np.mean(fC02_lon), np.mean(fC02_lat)]))

np.save('l4_mod_temp.npy', mod_temp[:,l4_mean])
np.save('l4_mod_sal.npy', mod_sal[:,l4_mean])
np.save('l4_mod_dt.npy', mod_dt)

##############################################################
# Get the reference salinity at each observation

ref_sal = []

for i, this_dt in enumerate(fC02_dt):
    mod_time_ind = np.argmin(np.asarray([np.abs((this_mod_dt - this_dt).total_seconds())
									for this_mod_dt in mod_dt]))
    mod_node_ind = tamar_grid.closest_node([fC02_lon[i], fC02_lat[i]])

    ref_sal.append(mod_sal[mod_time_ind, mod_node_ind])

ref_sal = np.squeeze(np.asarray(ref_sal))


###############################################################
# Extend the obs series to every model timestep - try linear interpolation first
ref_date = dt.datetime(2016,1,1)
mod_time_s = np.asarray([(this_t - ref_date).total_seconds() for this_t in mod_dt])
obs_time_s = np.asarray([(this_t - ref_date).total_seconds() for this_t in fC02_dt])

interp_ref_sal = np.interp(mod_time_s, obs_time_s, ref_sal)
interp_pC02 = np.interp(mod_time_s, obs_time_s, pC02) 
interp_pC02_atmos = np.interp(mod_time_s, obs_time_s, pC02_atmos)

"""
fig, ax = plt.subplots(figsize=[12,9])
ax2 = ax.twinx()

ax.plot(fC02_dt, fC02, c='C0')
ax.scatter(fC02_dt, fC02, c='C0')

ax2.scatter(fC02_dt, ref_sal, c='darkorange')
ax2.plot(mod_dt, mod_sal[:,l4_mean], c='bisque')
ax2.plot(mod_dt, interp_ref_sal, c='red')
ax2.plot(buoy_dt_ts, buoy_sal, c='lightgray')
"""

delta_sal = np.tile(interp_ref_sal[:,np.newaxis], [1,mod_sal.shape[1]]) - mod_sal
pC02_node = np.tile(interp_pC02[:,np.newaxis], [1,mod_sal.shape[1]]) - 39.83*delta_sal
deltaC02_node = np.tile(interp_pC02_atmos[:,np.newaxis], [1,mod_sal.shape[1]]) - pC02_node


###############################################################
# Get the region of validity

sal_valid = np.logical_and(delta_sal > salinity_below_lim, delta_sal <= salinity_above_lim)

if static_depth:
	dep_valid = np.tile((tamar_grid.grid.h >= bathy_limit)[np.newaxis,:], [sal_valid.shape[0], 1])
else:
	dep_valid = tamar_grid.grid.h  + mod_zeta >= bathy_limit

valid_depsal = np.logical_and(sal_valid, dep_valid)
valid = np.logical_and(box_valid, valid_depsal)

#all_valid = np.all(valid, axis=0)
val_perc = np.sum(valid, axis=0)/valid.shape[0]
all_valid = val_perc >= time_in_lim


# Get integrated flux under different wind conditions
grid_flux, grid_fmo, grid_fmolhr = fcf.flux_c02_calc(cnst_wind, mod_temp, deltaC02_node, mod_sal) 


tot_flux = grid_fmolhr*tamar_grid.grid.art1
flux_int = np.sum(tot_flux[:,all_valid], axis=1)


# Wind spd
buoy_time_s = np.asarray([(this_t - ref_date).total_seconds() for this_t in buoy_dt])
buoy_wind_interp = np.tile(np.interp(mod_time_s, buoy_time_s, wind_spd)[:,np.newaxis], [1,deltaC02_node.shape[1]])
grid_flux_var, grid_fmo_var, grid_fmolhr_var = fcf.flux_c02_calc(buoy_wind_interp, mod_temp, deltaC02_node, mod_sal)

tot_flux_var = grid_fmolhr_var*tamar_grid.grid.art1
flux_int_var = np.sum(tot_flux_var[:,all_valid], axis=1)

choose_out = np.arange(0,5881,24)

out_data = {'grid_flux':grid_flux[choose_out,:], 'time_dt':mod_dt, 'subset':choose_out, 'all_valid':all_valid, 'flux_int':flux_int}
out_data_var = {'grid_flux':grid_flux_var[choose_out,:], 'time_dt':mod_dt, 'subset':choose_out, 'all_valid':all_valid, 'flux_int':flux_int_var}

np.save('out_data_{}_7_1_fmolhr.npy'.format(area_ident), out_data)
np.save('out_data_{}_var_fmolhr.npy'.format(area_ident), out_data_var)
