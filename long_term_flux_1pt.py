import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

import flux_calc_funcs as fcf

wind_spd_cst = 7.1
area_file = 'out_data_all_var_fmolhr.npy'

##################################################################
# Get the buoy wind, temp, sal
wrf_l4_data = np.load('wrf_l4_wnd.npy', allow_pickle=True)
buoy_dt = wrf_l4_data[0]
buoy_dt, buoy_ind = np.unique(buoy_dt, return_index=True)

wind_spd = np.asarray(wrf_l4_data[1], dtype=float)[buoy_ind]



###############################################################
# construct interpolated PC02 etc from Rich

l4_fC02_file = 'Quest_2016_fCO2_L4_only.csv'

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
# Get landschutz data
ls_file = 'MPI-ULB-SOM_FFN_clim.nc'
ls_nc = nc.Dataset(ls_file)

bounds = [[-4.2, -4.01], [50.1, 50.2]]

ls_lon = ls_nc.variables['lon'][:]
ls_lat = ls_nc.variables['lat'][:]

choose_lon = np.logical_and(ls_lon >= bounds[0][0], ls_lon <= bounds[0][1])
choose_lat = np.logical_and(ls_lat >= bounds[1][0], ls_lat <= bounds[1][1])

landsch_pC02 = np.squeeze(ls_nc['pco2'][:,choose_lat, choose_lon])

landsch_dt = np.asarray([dt.datetime(2016,m,15) for m in np.arange(1,13)])


##############################################################
# Get mod l4 sal and temp
mod_l4_sal = np.load('l4_mod_sal.npy')
mod_l4_temp = np.load('l4_mod_temp.npy')
mod_dt = np.load('l4_mod_dt.npy', allow_pickle=True)


###############################################################
# Extend the obs series to hourly timesteps

ref_date = dt.datetime(2016,1,1)
mod_time_s = np.asarray([(this_t - ref_date).total_seconds() for this_t in mod_dt])
obs_time_s = np.asarray([(this_t - ref_date).total_seconds() for this_t in fC02_dt])
landsch_time_s = np.asarray([(this_t - ref_date).total_seconds() for this_t in landsch_dt])

interp_pC02 = np.interp(mod_time_s, obs_time_s, pC02) 
interp_pC02_atmos = np.interp(mod_time_s, obs_time_s, pC02_atmos)
interp_landsch_pC02 =  np.interp(mod_time_s, landsch_time_s, landsch_pC02)

obs_dc02 = interp_pC02_atmos - interp_pC02
landsch_dc02 = interp_pC02_atmos - interp_landsch_pC02

buoy_time_s = np.asarray([(this_t - ref_date).total_seconds() for this_t in buoy_dt])
buoy_wind_interp = np.interp(mod_time_s, buoy_time_s, wind_spd)


###############################################################
# Get integrated flux under different wind conditions
_, _, obs_flux_fmohl_var = fcf.flux_c02_calc(buoy_wind_interp, mod_l4_temp, obs_dc02, mod_l4_sal)
_, _, obs_flux_fmohl_cst = fcf.flux_c02_calc(wind_spd_cst, mod_l4_temp, obs_dc02, mod_l4_sal)
_, _, landsch_flux_fmohl_var = fcf.flux_c02_calc(buoy_wind_interp, mod_l4_temp, landsch_dc02, mod_l4_sal)
_, _, landsch_flux_fmohl_cst = fcf.flux_c02_calc(wind_spd_cst, mod_l4_temp, landsch_dc02, mod_l4_sal)

# Get grid area to multiply up
area_data = np.load(area_file, allow_pickle=True).item()
valid_area = np.sum(area_data['grid_area'][area_data['all_valid']])

out_data = {'obs_var_flux':obs_flux_fmohl_var*valid_area, 'obs_cst_flux':obs_flux_fmohl_cst*valid_area, 'landsch_var_flux':landsch_flux_fmohl_var*valid_area, 'landsch_cst_flux':landsch_flux_fmohl_cst*valid_area, 'time_dt':mod_dt}

np.save('pt_flux_fmolh.npy', out_data)

