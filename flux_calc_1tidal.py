import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

import PyFVCOM as pf

import flux_calc_funcs as fcf

months = [6,7,8,9]
len_out = 14

for month in months:
#    mod_basedir = '/data/sthenno1/backup/mbe/Models/FVCOM/tamar_v2/output/phys_only_2017/2016'
    #mod_basedir = '/data/sthenno1/scratch/mbe/Models_scratch/FVCOM/tamar_v2/phys_2017_sal35_10/2016'
    mod_basedir = '/data/proteus1/scratch/mbe/Models_scratch/FVCOM/tamar_v2/bias_corr2/2016'

    mod_filestr = '{}/{:02d}/tamar_v2_0001.nc'.format(mod_basedir, month)

    #mod_fr = pf.read.FileReader(mod_filestr, ['temp', 'salinity', 'uwind_speed', 'vwind_speed'])

    mod_fr =  pf.read.FileReader(mod_filestr, ['temp', 'salinity', 'zeta'])

    transfer_dir = '/users/modellers/mbe/transfer/'

    data_dir = '/data/sthenno1/backup/mbe/Experiments/Rich_paper/Data/'
    l4_fC02_file = '{}Quest_2016_fCO2_L4_only.csv'.format(data_dir)

    # Get the C02 measurements

    fC02_data = np.loadtxt(l4_fC02_file, delimiter=',', skiprows=1)
    fC02_dt_bst = np.asarray([dt.datetime(int(this_row[0]), int(this_row[1]), int(this_row[2]), int(this_row[3]), int(this_row[4]), int(this_row[5])) for this_row in fC02_data])
    fC02_dt = np.asarray([this_dt - dt.timedelta(hours=1) for this_dt in fC02_dt_bst])

    fC02_day = np.unique(np.asarray([dt.datetime(2016, this_dt.month, this_dt.day) for this_dt in fC02_dt]))

    fC02_month = np.asarray([this_dt.month for this_dt in fC02_day])
    chosen_dates = fC02_day[fC02_month == month]

    """
    # Get buoy data for wind speed
    this_year = 2016
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

    buoy_dt = buoy_dt[~remove_duff]
    wind_spd = np.sqrt(wind_u[~remove_duff]**2 + wind_v[~remove_duff]**2)
    """
    ### Use WRF wind
    wrf_l4_data = np.load('wrf_l4_wnd.npy', allow_pickle=True)
    buoy_dt = wrf_l4_data[0] 
    wind_spd = np.asarray(wrf_l4_data[1], dtype=float)

    ### Define the exlusions
    bathy_limit = 10
    salinity_below_lim = -4
    salinity_above_lim = 0.5


    # Run for each experiment date
    transfer_dir = '/users/modellers/mbe/transfer/' # for plots


    for this_date in chosen_dates:
        choose_dt = np.logical_and(fC02_dt >= this_date, fC02_dt <= this_date + dt.timedelta(days=1))

        out_date_fmt = this_date.strftime('%Y-%m-%d')

        poss_dates = fC02_dt[choose_dt]
        mid_date = poss_dates[int(np.floor(len(poss_dates)/2))]

        fC02 = np.mean(fC02_data[choose_dt,8])
        pC02 = np.mean(fC02_data[choose_dt,6])
        pC02_atmos = np.mean(fC02_data[choose_dt,7])

        fC02_lon = np.mean(fC02_data[choose_dt,13])
        fC02_lat = np.mean(fC02_data[choose_dt,12])
        l4_ll = [fC02_lon, fC02_lat]

        mod_l4 = mod_fr.closest_node(l4_ll)[0]
        mod_bathy = mod_fr.grid.h

        mod_times = np.asarray([mid_date + dt.timedelta(hours=int(this_h)) for this_h in np.arange(0,len_out,1)])

        mod_zeta = []
        out_ref_sal = []
        out_pc02 = []
        out_flux = []
        out_dc02 = []
        out_fmo = []
        out_fmolhr = []
        out_delta_sal = []
        out_temp = []
        out_wnd_spd = []
        out_flux_fixed = []
        out_flux_fmohl_fixed = []

        all_wind_spd = []

        for this_time in mod_times:
            closest_buoy_dt = np.argmin([(this_time - buoy_dt).total_seconds()**2 for buoy_dt in buoy_dt])
            all_wind_spd.append(wind_spd[closest_buoy_dt])
        
        mean_wnd = np.nanmean(np.asarray(all_wind_spd))

        first_time = True

        for this_time in mod_times:    

            mod_calc_time = mod_fr.closest_time(this_time)

            if first_time:
                ref_sal = mod_fr.data.salinity[mod_calc_time, 0, mod_l4]
                first_time = False

            delta_sal = ref_sal - mod_fr.data.salinity[mod_calc_time, 0, :]
        
            mod_zeta.append(mod_fr.data.zeta[mod_calc_time,:])

            pC02_node = pC02 - 39.83*delta_sal
            deltaC02_node = pC02_atmos - pC02_node

            closest_buoy_dt = np.argmin([(this_time - buoy_dt).total_seconds()**2 for buoy_dt in buoy_dt])

            # Run the flux calc
            grid_flux, grid_fmo, grid_fmolhr = fcf.get_mod_flux(mod_fr, this_time, deltaC02_node, wind_speed=wind_spd[closest_buoy_dt])

            time_ind = mod_fr.closest_time(this_time)
            out_ref_sal.append(ref_sal)
            out_pc02.append(pC02_node)
            out_flux.append(grid_flux)
            out_dc02.append(deltaC02_node)
            out_fmo.append(grid_fmo)
            out_fmolhr.append(grid_fmolhr)
            out_delta_sal.append(delta_sal)
            out_temp.append(mod_fr.data.temp[time_ind,0,:])
            out_wnd_spd.append(wind_spd[closest_buoy_dt])
            
            fixed_grid_flux, fixed_grid_fmo, fixed_grid_fmolhr = fcf.get_mod_flux(mod_fr, this_time, deltaC02_node, wind_speed=mean_wnd)
            out_flux_fixed.append(fixed_grid_flux)
            out_flux_fmohl_fixed.append(fixed_grid_fmolhr)

        save_list = {'date_str':out_date_fmt, 'ref_sal':np.asarray(out_ref_sal), 'dsal':np.asarray(out_delta_sal),'pc02':np.asarray(out_pc02), 'dc02':np.asarray(out_dc02), 'flux':np.asarray(out_flux), 'fmo':grid_fmo, 'fmohl':np.asarray(out_fmolhr), 'mod_zeta':np.asarray(mod_zeta), 'mod_h':mod_fr.grid.h, 'l4_loc':l4_ll, 'mod_times':mod_times, 'mod_temp':out_temp, 'wind_speed':out_wnd_spd, 'mean_wnd':mean_wnd, 'flux_fixed':out_flux_fixed, 'flux_fmohl_fixed':out_flux_fmohl_fixed}

        np.save('tidal_flux_{}.npy'.format(out_date_fmt), save_list, allow_pickle=True)

