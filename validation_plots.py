import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import glob as gb
from cmocean import cm

import PyFVCOM as pf
import matplotlib.gridspec as gridspec


exp_dates = {'low':dt.datetime(2016,7,7), 'flood':dt.datetime(2016,6,15), 'high':dt.datetime(2016,6,30), 'ebb':dt.datetime(2016,6,10)}
exp_states = {dt.datetime(2016,7,7):'low', dt.datetime(2016,6,15):'flood', dt.datetime(2016,6,30):'high', dt.datetime(2016,6,10):'ebb'}

tamar_fstr = 'out.nc'
tamar_fr = pf.read.FileReader(tamar_fstr)

l4_ll = [-4.13, 50.15]
l4_node = tamar_fr.closest_node(l4_ll)[0]

rect_area = [[-4.25, 50.12] , [-4.05,50.5]]

min_dsal = -0.1
max_dsal = 4
depth_lim = 10

spatial_sal_cmap = cm.haline

flux_int_cmap = 'viridis' 
flux_int_sal_subplt = False

save_pdf = True

########################################################################################
# Plot buoy and quest data

def closest_times(base_dt, match_dt):
    match_inds = []
    for this_dt in match_dt:
        dist_s = np.asarray([np.abs((this_dt - td).total_seconds()) for td in base_dt])
        match_inds.append(np.argmin(dist_s))
    return np.asarray(match_inds)

buoy_dt, buoy_temp, buoy_sal= np.load('buoy_temp.npy', allow_pickle=True)

quest_dict = np.load('all_quest_match.npy', allow_pickle=True).item()

all_mod_delta_s = []
all_obs_delta_s = []

new_plot_dict = {}

for out_date_fmt, quest_data in quest_dict.items():
    print(out_date_fmt)
    quest_dt = quest_data[0]
    if quest_dt >= dt.datetime(2016,6,1) and quest_dt<dt.datetime(2016,8,1):
        mod_dt = quest_data[1]
        mod_ll = quest_data[2]
        mod_t_series = quest_data[3]
        mod_s_series = quest_data[4]
        mod_l4_series = quest_data[5]
        obs_t = quest_data[6]
        obs_s = quest_data[7]
        obs_dt = quest_data[8]
        obs_ll = quest_data[9]

        choose_dt = np.logical_and(obs_dt > quest_dt, obs_dt < quest_dt + dt.timedelta(days=1))

        obs_s = obs_s[choose_dt]
        obs_dt = obs_dt[choose_dt]
        obs_ll = obs_ll[choose_dt]

        buoy_choose = np.logical_and(buoy_dt <= np.max(obs_dt) + dt.timedelta(days=1), buoy_dt >= np.min(obs_dt) - dt.timedelta(days=1))
        buoy_choose_dt = buoy_dt[buoy_choose]
        buoy_choose_sal = buoy_sal[buoy_choose]

        l4_buoy_dts = closest_times(buoy_choose_dt, obs_dt)
        l4_s = buoy_choose_sal[l4_buoy_dts]

        mod_delta_s = mod_l4_series - mod_s_series
        obs_delta_s = l4_s - obs_s

        all_mod_delta_s.append(mod_delta_s)
        all_obs_delta_s.append(obs_delta_s)

        x_lims = [-4.35, -4.05]
        y_lims = [50.2, 50.42]

        if quest_dt in list(exp_dates.values()):
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[16,10])
            fig.suptitle('Salinity comp')

            ax1.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray', linewidth=0.5)
            ax1.scatter(obs_ll[:,0], obs_ll[:,1], c='r', zorder=2)
            ax1.set_xlim(x_lims)
            ax1.set_ylim(y_lims)

            ax2.plot(obs_dt, obs_delta_s, c='r', linewidth=2)
            ax2.plot(obs_dt, mod_delta_s, c='b', linewidth=2)

            plt.tight_layout()
            plt.savefig('sal_comp_2{}.png'.format(out_date_fmt), dpi=180)
            if save_pdf:
                 plt.savefig('sal_comp_2{}.pdf'.format(out_date_fmt), dpi=60)
            plt.close()

            new_plot_dict[exp_states[quest_dt]] = [obs_dt, obs_delta_s, mod_delta_s]

fig3 = plt.figure(figsize=(16,20), constrained_layout=True)
gs = fig3.add_gridspec(2, 2)

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['low']
f3_ax1 = fig3.add_subplot(gs[0, 0])
f3_ax1.plot(obs_dt, obs_delta_s, c='r', linewidth=2)
f3_ax1.plot(obs_dt, mod_delta_s, c='b', linewidth=2)

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['flood']
f3_ax2 = fig3.add_subplot(gs[0, 1])
f3_ax2.plot(obs_dt, obs_delta_s, c='r', linewidth=2)
f3_ax2.plot(obs_dt, mod_delta_s, c='b', linewidth=2)

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['high']
f3_ax3 = fig3.add_subplot(gs[1, 0])
f3_ax3.plot(obs_dt, obs_delta_s, c='r', linewidth=2)
f3_ax3.plot(obs_dt, mod_delta_s, c='b', linewidth=2)

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['ebb']
f3_ax4 = fig3.add_subplot(gs[1, 1])
f3_ax4.plot(obs_dt, obs_delta_s, c='r', linewidth=2)
f3_ax4.plot(obs_dt, mod_delta_s, c='b', linewidth=2)

plt.tight_layout()
plt.savefig('quest_salinity_line.png', dpi=180)
if save_pdf:
    plt.savefig('quest_salinity_line.pdf', dpi=60)
plt.close()


all_mod_delta_s = np.hstack(all_mod_delta_s)
all_obs_delta_s = np.hstack(all_obs_delta_s)

# Scatter plot of all quest data
"""
rem_lim = [-0.1, 1.5]
rem_data1 = np.logical_or(all_mod_delta_s > rem_lim[1], all_obs_delta_s > rem_lim[1])
rem_data2 = np.logical_or(all_mod_delta_s < rem_lim[0], all_obs_delta_s < rem_lim[0])
rem_data = np.logical_or(rem_data1, rem_data2)

plt.hist2d(all_mod_delta_s[~rem_data], all_obs_delta_s[~rem_data],bins=100, cmap='Blues', vmax=50)
plt.plot(rem_lim, rem_lim, c='orange')
"""

plt_lims = [0,2.5]
plt.figure(figsize=[10,12])
plt.scatter(all_mod_delta_s, all_obs_delta_s)
plt.xlim([0,2.5])
plt.ylim([0,2.5])
plt.plot([0,2.5], [0,2.5],c='orange')
plt.tight_layout()
plt.xlabel('Model delta salinity')
plt.ylabel('Underway delta salinity')
plt.tight_layout()
plt.savefig('quest_salinity_scatter.png', dpi=180)
if save_pdf:
    plt.savefig('quest_salinity_scatter.pdf', dpi=60)

plt.close()

# L4 CTD plot



#######################################################################################
# Plot spatial varying data for each exp
sal_min = 32
sal_max = 35.5

for this_date in exp_dates.values():
    this_flux_file = 'tidal_flux_{}.npy'.format(this_date.strftime('%Y-%m-%d'))

    flux_dict = np.load(this_flux_file, allow_pickle=True).item()

    sal = flux_dict['ref_sal'][0] - flux_dict['dsal'][0,:]
    dsal = flux_dict['dsal'][0,:]
    pc02 = flux_dict['dc02'][0,:]
    flux = np.asarray(flux_dict['flux_fixed'])[0,:]

    sal_valid = np.logical_and(flux_dict['dsal'][0,:] > min_dsal, flux_dict['dsal'][0,:] <= max_dsal)
    dep_valid = flux_dict['mod_h'] >= depth_lim

    box_valid = np.logical_and(np.logical_and(tamar_fr.grid.lon >= rect_area[0][0], tamar_fr.grid.lon <= rect_area[1][0]), np.logical_and(tamar_fr.grid.lat >= rect_area[0][1], tamar_fr.grid.lat <= rect_area[1][1]))

    valid = np.logical_and(sal_valid, dep_valid)
    valid = np.logical_and(box_valid, valid)

    all_valid_nodes = np.where(valid)[0]
    all_valid_ele = np.where(np.isin(tamar_fr.grid.triangles,all_valid_nodes))[0]

    all_valid_c = np.zeros(len(tamar_fr.grid.lonc), dtype=bool)
    all_valid_c[all_valid_ele] = True


    fig3 = plt.figure(figsize=(16,20), constrained_layout=True)
    gs = fig3.add_gridspec(2, 2)

    f3_ax1 = fig3.add_subplot(gs[0, 0])    
    f3_ax1.set_facecolor('olive')
    f3_ax1.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    f3_ax1.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(sal.shape))
    f3_col = f3_ax1.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, sal, zorder=2, cmap=cm.haline, vmin=sal_min, vmax=sal_max)
    f3_ax1.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax1.set_xlim([-4.28, -4.02])
    f3_ax1.set_ylim([50.07, 50.4])

    cax1 = f3_ax1.inset_axes([0.1, 0.06, 0.8, 0.03])
    fig3.colorbar(f3_col, cax=cax1, orientation='horizontal')

    f3_ax2 = fig3.add_subplot(gs[0, 1])
    f3_ax2.set_facecolor('olive')
    f3_ax2.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    f3_ax2.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(sal.shape), cmap='bwr')
    f3_col2 = f3_ax2.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, dsal, zorder=2, cmap=cm.haline)
    f3_ax2.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax2.set_xlim([-4.28, -4.02])
    f3_ax2.set_ylim([50.07, 50.4])

    cax2 = f3_ax2.inset_axes([0.1, 0.06, 0.8, 0.03])
    fig3.colorbar(f3_col2, cax=cax2, orientation='horizontal')


    f3_ax3 = fig3.add_subplot(gs[1, 0])
    f3_ax3.set_facecolor('olive')
    f3_ax3.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    f3_ax3.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(sal.shape), vmin=-1, vmax=1, cmap='bwr')
    f3_col3 = f3_ax3.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, pc02, mask=~all_valid_c, zorder=2, cmap=flux_int_cmap)
    f3_ax3.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax3.set_xlim([-4.28, -4.02])
    f3_ax3.set_ylim([50.07, 50.4])

    cax3 = f3_ax3.inset_axes([0.1, 0.06, 0.8, 0.03])
    fig3.colorbar(f3_col3, cax=cax3, orientation='horizontal')


    f3_ax3 = fig3.add_subplot(gs[1, 1])
    f3_ax3.set_facecolor('olive')
    f3_ax3.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    f3_ax3.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(sal.shape), vmin=-1, vmax=1, cmap='bwr')
    f3_col3 = f3_ax3.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, flux, mask=~all_valid_c, zorder=2, cmap=flux_int_cmap)
    f3_ax3.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax3.set_xlim([-4.28, -4.02])
    f3_ax3.set_ylim([50.07, 50.4])

    cax3 = f3_ax3.inset_axes([0.1, 0.06, 0.8, 0.03])
    fig3.colorbar(f3_col3, cax=cax3, orientation='horizontal')

    plt.tight_layout()
    fig3.savefig('dc02_{}.png'.format(this_date.strftime('%Y-%m-%d')), dpi=180)
    if save_pdf:
        fig3.savefig('dc02_{}.pdf'.format(this_date.strftime('%Y-%m-%d')), dpi=60)
    plt.close() 

#######################################################################################
# Plot the four tidal cycle 




tstp1 = 0
tstp2 = 3
tstp3 = 6
tstp4 = 9

end_stp = 14


for this_date in exp_dates.values():
    this_flux_file = 'tidal_flux_{}.npy'.format(this_date.strftime('%Y-%m-%d'))

    flux_dict = np.load(this_flux_file, allow_pickle=True).item()

    sal_valid = np.logical_and(flux_dict['dsal'] > min_dsal, flux_dict['dsal'] <= max_dsal)
    dep_valid = flux_dict['mod_h'] >= depth_lim
    zeta_valid = flux_dict['mod_h'] + flux_dict['mod_zeta'] >= depth_lim

    box_valid = np.logical_and(np.logical_and(tamar_fr.grid.lon >= rect_area[0][0], tamar_fr.grid.lon <= rect_area[1][0]), np.logical_and(tamar_fr.grid.lat >= rect_area[0][1], tamar_fr.grid.lat <= rect_area[1][1]))
    box_valid = np.tile(box_valid[np.newaxis,:], [zeta_valid.shape[0],1])

    valid = np.logical_and(sal_valid, np.tile(dep_valid, [sal_valid.shape[0], 1]))
    valid = np.logical_and(box_valid, valid)
    all_valid = np.all(valid, axis=0)


    all_valid_nodes = np.where(all_valid)[0]
    all_valid_ele = np.where(np.isin(tamar_fr.grid.triangles,all_valid_nodes))[0]

    all_valid_c = np.zeros(len(tamar_fr.grid.lonc), dtype=bool)
    all_valid_c[all_valid_ele] = True

#    flux_rate = flux_dict['flux']
#    tot_flux = flux_dict['flux']*tamar_fr.grid.art1
#    flux_int = np.sum(tot_flux[:,all_valid], axis=1)

    flux_rate_cnst = np.asarray(flux_dict['flux_fixed'])
    tot_flux_cnst = flux_rate_cnst*tamar_fr.grid.art1
    flux_int_cnst = np.sum(tot_flux_cnst[:,all_valid], axis=1)


    flux_max = np.max(flux_rate_cnst[:,all_valid])
    flux_min = np.min(flux_rate_cnst[:,all_valid])

    sal = np.tile(flux_dict['ref_sal'][:,np.newaxis], [1,len(all_valid)]) - flux_dict['dsal']
    sal_min = 34.5
    sal_max = np.max(sal)

    fig3 = plt.figure(figsize=(14,24), constrained_layout=True)
    gs = fig3.add_gridspec(9, 2)

    f3_ax1 = fig3.add_subplot(gs[0:2, :])
    #f3_ax1.plot(flux_dict['mod_times'][0:end_stp], flux_int[0:end_stp], c='crimson')
    f3_ax1.plot(flux_dict['mod_times'][0:end_stp], flux_int_cnst[0:end_stp],c='darkblue')

    f3_ax1.scatter(flux_dict['mod_times'][tstp1], flux_int_cnst[tstp1], c='r')
    f3_ax1.scatter(flux_dict['mod_times'][tstp2], flux_int_cnst[tstp2], c='r')
    f3_ax1.scatter(flux_dict['mod_times'][tstp3], flux_int_cnst[tstp3], c='r')
    f3_ax1.scatter(flux_dict['mod_times'][tstp4], flux_int_cnst[tstp4], c='r')

    f3_ax2 = fig3.add_subplot(gs[2, :])
    f3_ax2_wnd = f3_ax2.twinx()

    f3_ax2.plot(flux_dict['mod_times'][0:end_stp], flux_dict['mod_zeta'][0:end_stp,l4_node], c='darkgrey')

    #f3_ax2_wnd.plot(flux_dict['mod_times'][0:end_stp], flux_dict['wind_speed'][0:end_stp], c='crimson')
    f3_ax2_wnd.plot(flux_dict['mod_times'][0:end_stp], flux_dict['mean_wnd'] * np.ones(end_stp), c='darkblue', linestyle='--')
    f3_ax2.scatter(flux_dict['mod_times'][tstp1], flux_dict['mod_zeta'][tstp1,l4_node], c='k', zorder=2)
    f3_ax2.scatter(flux_dict['mod_times'][tstp2], flux_dict['mod_zeta'][tstp2,l4_node], c='k', zorder=2)
    f3_ax2.scatter(flux_dict['mod_times'][tstp3], flux_dict['mod_zeta'][tstp3,l4_node], c='k', zorder=2)
    f3_ax2.scatter(flux_dict['mod_times'][tstp4], flux_dict['mod_zeta'][tstp4,l4_node], c='k', zorder=2)

    

    f3_ax3 = fig3.add_subplot(gs[3:6, 0])
    f3_ax3.set_facecolor('olive')
    f3_ax3.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    f3_ax3.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(tot_flux_cnst[tstp1,:].shape), vmin=-1, vmax=1, cmap='bwr')
    f3_ax3.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, flux_rate_cnst[tstp1,:], mask=~all_valid_c, zorder=2, vmin=flux_min, vmax=flux_max, cmap=flux_int_cmap)
    f3_ax3.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax3.set_xlim([-4.28, -4.02])
    f3_ax3.set_ylim([50.07, 50.4])

    f3_ax4 = fig3.add_subplot(gs[3:6, 1])
    f3_ax4.set_facecolor('olive')
    f3_ax4.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    f3_ax4.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(tot_flux_cnst[tstp2,:].shape), vmin=-1, vmax=1, cmap='bwr')
    f3_ax4.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, flux_rate_cnst[tstp2,:], mask=~all_valid_c, zorder=2, vmin=flux_min, vmax=flux_max, cmap=flux_int_cmap)
    f3_ax4.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax4.set_xlim([-4.28, -4.02])
    f3_ax4.set_ylim([50.07, 50.4])

    f3_ax5 = fig3.add_subplot(gs[6:, 0])
    f3_ax5.set_facecolor('olive')
    f3_ax5.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    flux_trp = f3_ax5.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(tot_flux_cnst[tstp3,:].shape), vmin=-1, vmax=1, cmap='bwr')
    f3_ax5.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, flux_rate_cnst[tstp3,:], mask=~all_valid_c, zorder=2, vmin=flux_min, vmax=flux_max, cmap=flux_int_cmap)
    f3_ax5.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax5.set_xlim([-4.28, -4.02])
    f3_ax5.set_ylim([50.07, 50.4])

    f3_ax6 = fig3.add_subplot(gs[6:, 1])
    f3_ax6.set_facecolor('olive')
    f3_ax6.triplot(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, c='lightgray')
    f3_ax6.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, np.zeros(tot_flux_cnst[tstp3,:].shape), vmin=-1, vmax=1, cmap='bwr')
    flux_trp = f3_ax6.tripcolor(tamar_fr.grid.lon, tamar_fr.grid.lat, tamar_fr.grid.triangles, flux_rate_cnst[tstp3,:], mask=~all_valid_c, zorder=2, vmin=flux_min, vmax=flux_max, cmap=flux_int_cmap)
    f3_ax6.scatter(tamar_fr.grid.lon[l4_node], tamar_fr.grid.lat[l4_node], c='purple', zorder=2)
    f3_ax6.set_xlim([-4.28, -4.02])
    f3_ax6.set_ylim([50.07, 50.4])


    cax = f3_ax6.inset_axes([0.1, 0.06, 0.8, 0.03])
    fig3.colorbar(flux_trp, cax=cax, orientation='horizontal')

    plt.tight_layout()
    fig3.savefig('flux_tidal_{}.png'.format(flux_dict['date_str']), dpi=180)
    if save_pdf:
        fig3.savefig('flux_tidal_{}.pdf'.format(flux_dict['date_str']), dpi=60)
    plt.close()


