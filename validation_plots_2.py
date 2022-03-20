import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import glob as gb
import netCDF4 as nc
from cmocean import cm
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import pickle as pk
import PyFVCOM as pf


exp_dates = {'low':dt.datetime(2016,7,7), 'flood':dt.datetime(2016,6,15), 'high':dt.datetime(2016,6,30), 'ebb':dt.datetime(2016,6,10)}
exp_states = {dt.datetime(2016,7,7):'low', dt.datetime(2016,6,15):'flood', dt.datetime(2016,6,30):'high', dt.datetime(2016,6,10):'ebb'}

[tamar_fr_lon, tamar_fr_lonc, tamar_fr_lat, tamar_fr_latc, tamar_fr_triangles, tamar_fr_art1, l4_node] = np.load('fvcom_grid_data.npy', allow_pickle=True)

#mikes L4 coordinates were wrong redefine L4 node
l4_ll = [-4.221, 50.251]
tamar_fstr = 'out.nc'
tamar_fr = pf.read.FileReader(tamar_fstr)
l4_node = tamar_fr.closest_node(l4_ll)[0]


l4_label_pos = [-4.211, 50.251]
l4_node_label_pos = tamar_fr.closest_node(l4_label_pos)[0]

#area label
text_label_pos = [-4.165, 50.25]
text_node_label_pos = tamar_fr.closest_node(text_label_pos)[0]

rect_area = [[-4.25, 50.32] , [-4.05,50.50]]

#this is the resticted area north of Penlee
# rect_area = [[-4.193, 50.319] , [-4.05,50.50]]


max_dsal = +0.1
min_dsal = -4
depth_lim = 10

# spatial_sal_cmap = cm.haline

flux_int_cmap = 'viridis' 
flux_int_sal_subplt = False

save_pdf = False

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
            """
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=[16,10])
            fig.suptitle('Salinity comp')

            ax1.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray', linewidth=0.5)
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
            """
            new_plot_dict[exp_states[quest_dt]] = [obs_dt, obs_delta_s, mod_delta_s]


HM_Fmt = mdates.DateFormatter('%H:%M')
plt.rcParams.update({'font.size': 22})
fig3 = plt.figure(figsize=(22,20), constrained_layout=True)
gs = fig3.add_gridspec(2, 2)

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['low']
f3_ax1 = fig3.add_subplot(gs[0, 0])
f3_ax1.plot(obs_dt, obs_delta_s,  linewidth=2.5,color='black')
f3_ax1.plot(obs_dt, mod_delta_s,  linewidth=2.5,color='orange')
f3_ax1.text(0.0, 1.03, '(a)',transform=f3_ax1.transAxes,fontsize=30)
f3_ax1.set_title('LW+ 0.5 hrs',fontsize=30)
f3_ax1.xaxis.set_major_formatter(HM_Fmt)
f3_ax1.tick_params(axis="x", labelsize=25) 
f3_ax1.tick_params(axis="y", labelsize=25) 
f3_ax1.set_xlabel(r'Time (HH:MM)',fontsize=30,color='black')
f3_ax1.set_ylabel(r'$\xi$S (PSU)',fontsize=30,color='black')
plt.rc('image', cmap='viridis')

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['flood']
f3_ax2 = fig3.add_subplot(gs[0, 1])
f3_ax2.plot(obs_dt, obs_delta_s,  linewidth=2.5,color='black')
f3_ax2.plot(obs_dt, mod_delta_s,  linewidth=2.5,color='orange')
f3_ax2.text(0.0, 1.03, '(b)',transform=f3_ax2.transAxes, fontsize=30)
f3_ax2.set_title('LW+ 3.5 hrs',fontsize=30)
f3_ax2.xaxis.set_major_formatter(HM_Fmt)
f3_ax2.tick_params(axis="x", labelsize=25) 
f3_ax2.tick_params(axis="y", labelsize=25)
f3_ax2.set_xlabel(r'Time (HH:MM)',fontsize=30,color='black')
f3_ax2.set_ylabel(r'$\xi$S (PSU)',fontsize=30,color='black')
plt.rc('image', cmap='viridis')

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['high']
f3_ax3 = fig3.add_subplot(gs[1, 0])
f3_ax3.plot(obs_dt, obs_delta_s,  linewidth=2.5,color='black')
f3_ax3.plot(obs_dt, mod_delta_s,  linewidth=2.5,color='orange')
f3_ax3.text(0.0, 1.03, '(c)',transform=f3_ax3.transAxes,fontsize=30)
f3_ax3.set_title('LW+ 6 hrs',fontsize=30)
f3_ax3.tick_params(axis="x", labelsize=25) 
f3_ax3.tick_params(axis="y", labelsize=25)
f3_ax3.xaxis.set_major_formatter(HM_Fmt)
f3_ax3.set_xlabel(r'Time (HH:MM)',fontsize=30,color='black')
f3_ax3.set_ylabel(r'$\xi$S (PSU)',fontsize=30,color='black')
plt.rc('image', cmap='viridis')

[obs_dt, obs_delta_s, mod_delta_s] = new_plot_dict['ebb']
f3_ax4 = fig3.add_subplot(gs[1, 1])
f3_ax4.plot(obs_dt, obs_delta_s,  linewidth=2.5,color='black')
f3_ax4.plot(obs_dt, mod_delta_s,  linewidth=2.5,color='orange')
f3_ax4.text(0.0, 1.03, '(d)',transform=f3_ax4.transAxes,fontsize=30)
f3_ax4.set_title('LW+ 8.75 hrs',fontsize=30)
f3_ax4.tick_params(axis="x", labelsize=25) 
f3_ax4.tick_params(axis="y", labelsize=25)
f3_ax4.xaxis.set_major_formatter(HM_Fmt)
f3_ax4.set_xlabel(r'Time (HH:MM)',fontsize=30,color='black')
f3_ax4.set_ylabel(r'$\xi$S (PSU)',fontsize=30,color='black')
plt.rc('image', cmap='viridis')


plt.tight_layout()
plt.savefig('quest_salinity_line.png', dpi=180)
if save_pdf:
    plt.savefig('quest_salinity_line.pdf', dpi=60)
plt.close()


all_mod_delta_s = np.hstack(all_mod_delta_s)
all_obs_delta_s = np.hstack(all_obs_delta_s)

import matplotlib.pyplot as plt
from scipy import stats
float_arr = np.vstack(all_obs_delta_s).astype(np.float32)
float_arr2 = float_arr.flatten()
#these are the full data!
slope, intercept, r_value, p_value, std_err = stats.linregress(all_mod_delta_s,float_arr2)
RMSE_SAL=np.sqrt(np.mean((float_arr2-all_mod_delta_s)**2))
R_squared_sal=r_value*r_value

# Scatter plot of all quest data

rem_lim = [-0.1, 4]
rem_data1 = np.logical_or(all_mod_delta_s > rem_lim[1], all_obs_delta_s > rem_lim[1])
rem_data2 = np.logical_or(all_mod_delta_s < rem_lim[0], all_obs_delta_s < rem_lim[0])
rem_data = np.logical_or(rem_data1, rem_data2)

plt.hist2d(all_mod_delta_s[~rem_data], all_obs_delta_s[~rem_data],bins=100, cmap='Blues', vmax=50)
plt.plot(rem_lim, rem_lim, c='black')


#VERIFIED that you get a lower RMSD using those salinity limits above.
RMSE_SAL_2=np.sqrt(np.mean((all_mod_delta_s[~rem_data]-all_obs_delta_s[~rem_data])**2))
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(all_mod_delta_s[~rem_data],float_arr2[~rem_data])
R_squared_sal2=r_value2*r_value2

plt_lims = [0,4]
plt.figure(figsize=[10,12])
plt.scatter(all_mod_delta_s, all_obs_delta_s)
plt.xlim([0,2.5])
plt.ylim([0,2.5])
plt.plot([0,2.5], [0,2.5],'--',c='black')

# Data for plotting
t = np.arange(0.0, 2.05, 0.01)
s = (t*slope)+intercept
s2 = (t*slope2)+intercept2
plt.plot(t,s2,'--',c='red')

plt.tight_layout()
plt.xlabel(r'Modelled $\xi$S (PSU)',fontsize=26,color='black')
plt.ylabel(r'Observed $\xi$S (PSU)',fontsize=26,color='black')
plt.tight_layout()
equation = 'y = ' + str(round(slope2,4)) + 'x' + str(round(intercept2,4))
plt.text(0.5, 2,equation,fontsize=26,color='red')
equation2 = 'R\u00b2 = '  + str(round(R_squared_sal2,4)) 
plt.text(0.5, 1.8,equation2,fontsize=26,color='red')
plt.savefig('quest_salinity_scatter.png', dpi=180)
if save_pdf:
    plt.savefig('quest_salinity_scatter.pdf', dpi=60)

plt.close()




# L4 CTD plot



#######################################################################################
# Plot spatial varying data for each exp
sal_min = 32.0
sal_max = 35.5

plt.rcParams.update({'font.size': 22})
for this_date in exp_dates.values():
    this_flux_file = 'tidal_flux_{}.npy'.format(this_date.strftime('%Y-%m-%d'))
    
    flux_dict = np.load(this_flux_file, allow_pickle=True).item()

    sal = flux_dict['ref_sal'][6] + flux_dict['dsal'][6,:]
    dsal = flux_dict['dsal'][6,:]
    pc02 = flux_dict['dc02'][6,:]
    fc02 = flux_dict['fc02'][6,:]
    # flux = -np.asarray(flux_dict['flux_fixed'])[0,:]# check what the units are of this!

    sal_valid = np.logical_and(flux_dict['dsal'][6,:] > min_dsal, flux_dict['dsal'][6,:] <= max_dsal)
    dep_valid = flux_dict['mod_h'] >= depth_lim

    box_valid = np.logical_and(np.logical_and(tamar_fr_lon >= rect_area[0][0], tamar_fr_lon <= rect_area[1][0]), np.logical_and(tamar_fr_lat >= rect_area[0][1], tamar_fr_lat <= rect_area[1][1]))

    valid = np.logical_and(sal_valid, dep_valid)
    valid = np.logical_and(box_valid, valid)

    all_valid_nodes = np.where(valid)[0]
    all_valid_ele = np.where(np.isin(tamar_fr_triangles,all_valid_nodes))[0]

    all_valid_c = np.zeros(len(tamar_fr_lonc), dtype=bool)
    all_valid_c[all_valid_ele] = True


    fig3 = plt.figure(figsize=(16,10), constrained_layout=True)
    gs = fig3.add_gridspec(1, 2)
    
    # fig3.subplots_adjust(top=0.98)
    # fig3.subplots_adjust(left=0.05)
    # fig3.subplots_adjust(right=0.98)
    # fig3.subplots_adjust(wspace=0.16)
    # fig3.subplots_adjust(hspace=0.2)
    # fig3.subplots_adjust(bottom=0.05)
    
    f3_ax1 = fig3.add_subplot(gs[0, 0]) 
    f3_ax1.set_facecolor('olive')
    f3_ax1.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    f3_ax1.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(sal.shape))
    f2_col = f3_ax1.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, sal, zorder=2, cmap=flux_int_cmap, vmin=sal_min, vmax=sal_max)
    # cmap='viridis_r'
    f3_ax1.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=150, c='white', zorder=2)
    f3_ax1_txt_l4 = f3_ax1.annotate( 'L4',[tamar_fr_lon[l4_node_label_pos], tamar_fr_lat[l4_node_label_pos]], c='white')
    f3_ax1_txt = f3_ax1.text(0.00, 1.03, '(a)', transform=f3_ax1.transAxes, fontsize='x-large')
    f3_ax1_txt.set_bbox(dict(facecolor='white', edgecolor='white'))
    f3_ax1.set_xlim([-4.28, -4.02])
    f3_ax1.set_ylim([50.23, 50.4])
    f3_ax1_xlab = f3_ax1.set_xlabel('$^{\circ}$ Lon')
    f3_ax1.xaxis.set_label_coords(-0.05, -0.02)
    f3_ax1_ylab = f3_ax1.set_ylabel('$^{\circ}$ Lat', rotation=0)
    f3_ax1.yaxis.set_label_coords(-0.08, 0.9)   
    
    cax1 = f3_ax1.inset_axes([0.1, -0.1, 0.8, 0.03])
    f3_cb1 = fig3.colorbar(f2_col, cax=cax1, orientation='horizontal')
    f3_cb1.set_label('Salinity (PSU)')
    
    f3_ax2 = fig3.add_subplot(gs[0, 1])
    f3_ax2.set_facecolor('olive')
    f3_ax2.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    f3_ax2.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(sal.shape), vmin=-1, vmax=1, cmap='bwr')
    f3_col2 = f3_ax2.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, fc02, mask=~all_valid_c, zorder=2, cmap=flux_int_cmap)
    f3_ax2.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=150, c='white', zorder=2)
    f3_ax2_txt_l4 = f3_ax2.annotate( 'L4',[tamar_fr_lon[l4_node_label_pos], tamar_fr_lat[l4_node_label_pos]], c='white')
    f3_ax2_txt = f3_ax2.text(0.00, 1.03, '(b)', transform=f3_ax2.transAxes, fontsize='x-large')
    f3_ax2_txt.set_bbox(dict(facecolor='white', edgecolor='white'))
    f3_ax2.set_xlim([-4.28, -4.02])
    f3_ax2.set_ylim([50.23, 50.4])
    f3_ax2_xlab = f3_ax2.set_xlabel('$^{\circ}$ Lon')
    f3_ax2.xaxis.set_label_coords(-0.05, -0.02)
    f3_ax2_ylab = f3_ax2.set_ylabel('$^{\circ}$ Lat', rotation=0)
    f3_ax2.yaxis.set_label_coords(-0.08, 0.9)  

    cax2 = f3_ax2.inset_axes([0.1, -0.1, 0.8, 0.03])
    f3_cb2 = fig3.colorbar(f3_col2, cax=cax2, orientation='horizontal')
    f3_cb2.set_label(r'$fCO_{2 (sw)}$' '(' '\u03BC' 'atm)')

    plt.tight_layout()
    fig3.savefig('dc02_{}.png'.format(this_date.strftime('%Y-%m-%d')), dpi=180)
    if save_pdf:
        fig3.savefig('dc02_{}.pdf'.format(this_date.strftime('%Y-%m-%d')), dpi=60)
    plt.close() 


    # mikes old plot - heavily edited above
    # fig3 = plt.figure(figsize=(16,20))
    # gs = fig3.add_gridspec(2, 2)

    # fig3.subplots_adjust(top=0.98)
    # fig3.subplots_adjust(left=0.05)
    # fig3.subplots_adjust(right=0.98)
    # fig3.subplots_adjust(wspace=0.16)
    # fig3.subplots_adjust(hspace=0.2)
    # fig3.subplots_adjust(bottom=0.05)

    # f3_ax1 = fig3.add_subplot(gs[0, 0])    
    # f3_ax1.set_facecolor('olive')
    # f3_ax1.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    # f3_ax1.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(sal.shape))
    # f3_col = f3_ax1.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, sal, zorder=2, cmap=cm.haline, vmin=sal_min, vmax=sal_max)
    # f3_ax1.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=150, c='red', zorder=2)
    # f3_ax1_txt_l4 = f3_ax1.text(0.48, 0.2, 'L4', transform=f3_ax1.transAxes, fontsize='x-large')
    # f3_ax1_txt = f3_ax1.text(0.07, 0.9, '(A)', transform=f3_ax1.transAxes, fontsize='x-large')
    # f3_ax1_txt.set_bbox(dict(facecolor='white', edgecolor='white'))
    # f3_ax1.set_xlim([-4.28, -4.02])
    # f3_ax1.set_ylim([50.07, 50.4])
    # f3_ax1_xlab = f3_ax1.set_xlabel('$^{\circ}$ Lon')
    # f3_ax1.xaxis.set_label_coords(-0.05, -0.02)
    # f3_ax1_ylab = f3_ax1.set_ylabel('$^{\circ}$ Lat', rotation=0)
    # f3_ax1.yaxis.set_label_coords(-0.08, 0.9)   
    
    # cax1 = f3_ax1.inset_axes([0.1, -0.1, 0.8, 0.03])
    # f3_cb1 = fig3.colorbar(f3_col, cax=cax1, orientation='horizontal')
    # f3_cb1.set_label('Salinity')

    # f3_ax2 = fig3.add_subplot(gs[0, 1])
    # f3_ax2.set_facecolor('olive')
    # f3_ax2.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    # f3_ax2.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(sal.shape), cmap='bwr')
    # f3_col2 = f3_ax2.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, dsal, zorder=2, cmap=cm.haline)
    # f3_ax2.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=150, c='red', zorder=2)
    # f3_ax2_txt = f3_ax2.text(0.07, 0.9, '(B)', transform=f3_ax2.transAxes, fontsize='x-large')
    # f3_ax2_txt.set_bbox(dict(facecolor='white', edgecolor='white'))
    

    # f3_ax2.set_xlim([-4.28, -4.02])
    # f3_ax2.set_ylim([50.07, 50.4])
    # f3_ax2_xlab = f3_ax2.set_xlabel('$^{\circ}$ Lon')
    # f3_ax2.xaxis.set_label_coords(-0.05, -0.02)
    # f3_ax2_ylab = f3_ax2.set_ylabel('$^{\circ}$ Lat', rotation=0)
    # f3_ax2.yaxis.set_label_coords(-0.08, 0.9)  

    # cax2 = f3_ax2.inset_axes([0.1, -0.1, 0.8, 0.03])
    # f3_cb2 = fig3.colorbar(f3_col2, cax=cax2, orientation='horizontal')
    # f3_cb2.set_label('$\Delta$ Salinity')
    

    # f3_ax3 = fig3.add_subplot(gs[1, 0])
    # f3_ax3.set_facecolor('olive')
    # f3_ax3.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    # f3_ax3.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(sal.shape), vmin=-1, vmax=1, cmap='bwr')
    # f3_col3 = f3_ax3.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, pc02, mask=~all_valid_c, zorder=2, cmap=flux_int_cmap)
    # f3_ax3.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=150, c='red', zorder=2)
    # f3_ax3_txt = f3_ax3.text(0.07, 0.9, '(C)', transform=f3_ax3.transAxes, fontsize='x-large')
    # f3_ax3_txt.set_bbox(dict(facecolor='white', edgecolor='white'))
    # f3_ax3.set_xlim([-4.28, -4.02])
    # f3_ax3.set_ylim([50.07, 50.4])
    # f3_ax3_xlab = f3_ax3.set_xlabel('$^{\circ}$ Lon')
    # f3_ax3.xaxis.set_label_coords(-0.05, -0.02)
    # f3_ax3_ylab = f3_ax3.set_ylabel('$^{\circ}$ Lat', rotation=0)
    # f3_ax3.yaxis.set_label_coords(-0.08, 0.9)  


    # cax3 = f3_ax3.inset_axes([0.1, -0.1, 0.8, 0.03])
    # f3_cb3 = fig3.colorbar(f3_col3, cax=cax3, orientation='horizontal')
    # f3_cb3.set_label('$\Delta$ fc02')
    

    # f3_ax4 = fig3.add_subplot(gs[1, 1])
    # f3_ax4.set_facecolor('olive')
    # f3_ax4.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    # f3_ax4.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(sal.shape), vmin=-1, vmax=1, cmap='bwr')
    # f3_col4 = f3_ax4.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, flux, mask=~all_valid_c, zorder=2, cmap=flux_int_cmap)
    # f3_ax4.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=150, c='red', zorder=2)
    # f3_ax4_txt = f3_ax4.text(0.07, 0.9, '(D)', transform=f3_ax4.transAxes, fontsize='x-large')
    # f3_ax4_txt.set_bbox(dict(facecolor='white', edgecolor='white'))
    # f3_ax4.set_xlim([-4.28, -4.02])
    # f3_ax4.set_ylim([50.07, 50.4])
    # f3_ax4_xlab = f3_ax4.set_xlabel('$^{\circ}$ Lon')
    # f3_ax4.xaxis.set_label_coords(-0.05, -0.02)
    # f3_ax4_ylab = f3_ax4.set_ylabel('$^{\circ}$ Lat', rotation=0)
    # f3_ax4.yaxis.set_label_coords(-0.08, 0.9)  

    # cax4 = f3_ax4.inset_axes([0.1, -0.1, 0.8, 0.03])
    # f3_cb4 = fig3.colorbar(f3_col4, cax=cax4, orientation='horizontal')
    # f3_cb4.set_label('C02 flux, mmol h$^{-1}$ m$^{2}$')
    

    # plt.tight_layout()
    # fig3.savefig('dc02_{}.png'.format(this_date.strftime('%Y-%m-%d')), dpi=180)
    # if save_pdf:
    #     fig3.savefig('dc02_{}.pdf'.format(this_date.strftime('%Y-%m-%d')), dpi=60)
    # plt.close() 

#######################################################################################
# Plot the four tidal cycle 

#stamps changes so is length 17 goes from - 3 hr to +14 hours
tstp1 = 0 # this is lw+ 0 hours
tstp2 = 3
tstp3 = 6
tstp4 = 9# this is lw +9 hours
end_stp = 12

plt.rcParams.update({'font.size': 22})
HM_Fmt = mdates.DateFormatter('%H:%M')

counter=0
fluxsave=dict()
meanwindsave=dict()
# Note this loop plots all four experiment dates but has been optimised for 2016-06-30
for this_date in exp_dates.values():
    this_flux_file = 'tidal_flux_{}.npy'.format(this_date.strftime('%Y-%m-%d'))
    counter=counter+1
    
    if counter==1:
       str_stp = 6 
       tstp1 = 6 # this is lw+ 0 hours
       tstp2 = 9
       tstp3 = 12
       tstp4 = 15# this is lw +9 hours
       end_stp = 17 
    if counter==2:
       str_stp = 3 
       tstp1 = 3 # this is lw+ 0 hours
       tstp2 = 6
       tstp3 = 9
       tstp4 = 12# this is lw +9 hours
       end_stp = 15 
    if counter==3:
       str_stp = 0  
       tstp1 = 0 # this is lw+ 0 hours
       tstp2 = 3
       tstp3 = 6
       tstp4 = 9# this is lw +9 hours
       end_stp = 12 
    if counter==4:
       str_stp = 6  
       tstp1 = 9 # this is lw+ 0 hours
       tstp2 = 12
       tstp3 = 15
       tstp4 = 6# this is lw +9 hours
       end_stp = 17 
    
    flux_dict = np.load(this_flux_file, allow_pickle=True).item()

    sal_valid = np.logical_and(flux_dict['dsal'] > min_dsal, flux_dict['dsal'] <= max_dsal)
    dep_valid = flux_dict['mod_h'] >= depth_lim
    zeta_valid = flux_dict['mod_h'] + flux_dict['mod_zeta'] >= depth_lim

    box_valid = np.logical_and(np.logical_and(tamar_fr_lon >= rect_area[0][0], tamar_fr_lon <= rect_area[1][0]), np.logical_and(tamar_fr_lat >= rect_area[0][1], tamar_fr_lat <= rect_area[1][1]))
    box_valid = np.tile(box_valid[np.newaxis,:], [zeta_valid.shape[0],1])

    valid = np.logical_and(sal_valid, np.tile(dep_valid, [sal_valid.shape[0], 1]))
    valid = np.logical_and(box_valid, valid)
    all_valid = np.all(valid, axis=0)


    all_valid_nodes = np.where(all_valid)[0]
    all_valid_ele = np.where(np.isin(tamar_fr_triangles,all_valid_nodes))[0]

    all_valid_c = np.zeros(len(tamar_fr_lonc), dtype=bool)
    all_valid_c[all_valid_ele] = True

#    flux_rate = flux_dict['flux']
#    tot_flux = flux_dict['flux']*tamar_fr.grid.art1
#    flux_int = np.sum(tot_flux[:,all_valid], axis=1)

    




    flux_rate_cnst = -np.asarray(flux_dict['flux_fmohl_fixed'])/1000# divide by 100 to get into m mol hr-1
    tot_flux_cnst = flux_rate_cnst*tamar_fr_art1
    flux_int_cnst = np.sum(tot_flux_cnst[:,all_valid], axis=1) #
    flux_int_cnst_mol = np.sum(tot_flux_cnst[:,all_valid], axis=1)/1000 #to get moles divide again


    vaid_area = round((np.sum(tamar_fr_art1[all_valid], axis=0))/1e6)#convert from m2 to km 2
    
    flux_rate_cnst_micromol = -np.asarray(flux_dict['flux_fmohl_fixed'])# in umol hr-1
    
    #this is just a stat neede for the paper
    if counter==1:
        fluxjuly7save=flux_int_cnst_mol[str_stp:end_stp]
    
    
    flux_max = np.max(flux_rate_cnst[:,all_valid])
    flux_min = np.min(flux_rate_cnst[:,all_valid])

    # save the integrated flux number at the three timestamps
    fluxsave[counter,0]=(flux_int_cnst_mol[tstp1])
    fluxsave[counter,1]=(flux_int_cnst_mol[tstp2])
    fluxsave[counter,2]=(flux_int_cnst_mol[tstp3])
    fluxsave[counter,3]=(flux_int_cnst_mol[tstp4])

    variable_win = flux_dict['mean_wnd']
    #meanwindsave[counter]=print(flux_dict['mean_wnd'])
    meanwindsave[counter]=variable_win


    fig3 = plt.figure(figsize=(16,24))
    gs = fig3.add_gridspec(46, 50)
    # fig3.subplots_adjust(top=0.98)
    # fig3.subplots_adjust(left=0.05)
    # fig3.subplots_adjust(right=0.98)
    # fig3.subplots_adjust(wspace=0.16)
    # fig3.subplots_adjust(hspace=0.2)
    # fig3.subplots_adjust(bottom=0.05)
    
    f3_ax1 = fig3.add_subplot(gs[0:4, 5:50])

    f3_ax1.plot(flux_dict['mod_times'][str_stp:end_stp], flux_dict['mod_zeta'][str_stp:end_stp,l4_node], c='black')
    
    f3_ax1.scatter(flux_dict['mod_times'][tstp1], flux_dict['mod_zeta'][tstp1,l4_node], s=150, c='r', zorder=2)
    f3_ax1.scatter(flux_dict['mod_times'][tstp2], flux_dict['mod_zeta'][tstp2,l4_node], s=150, c='r', zorder=2)
    f3_ax1.scatter(flux_dict['mod_times'][tstp3], flux_dict['mod_zeta'][tstp3,l4_node], s=150, c='r', zorder=2)
    f3_ax1.scatter(flux_dict['mod_times'][tstp4], flux_dict['mod_zeta'][tstp4,l4_node], s=150, c='r', zorder=2)

    f3_ax1_txt = f3_ax1.text(-0.25, 1.0, '(a)', transform=f3_ax1.transAxes, fontsize='x-large')
    f3_ax1.set_ylabel('Tidal \n height (m)')
    f3_ax1.xaxis.set_major_formatter(HM_Fmt)
    f3_ax1.set_xlabel('Time (HH:MM)')

    f3_ax2 = fig3.add_subplot(gs[6:10, 5:50])
    #f3_ax1.plot(flux_dict['mod_times'][0:end_stp], flux_int[0:end_stp], c='crimson')
    f3_ax2.plot(flux_dict['mod_times'][str_stp:end_stp], flux_int_cnst_mol[str_stp:end_stp],c='black')

    label_adj = [0.9*10**9, 0.7*10**9, 0.7*10**9, 0.6*10**9]
    label_x_adj = dt.timedelta(minutes=30)
    label_y_adj = 0.1*10**9
    label_yunder_adj = 0.12*10**9

    f3_ax2.scatter(flux_dict['mod_times'][tstp1], flux_int_cnst_mol[tstp1], s=150, c='r')
    f3_ax2.text(flux_dict['mod_times'][tstp1] + label_x_adj, flux_int_cnst_mol[tstp1], 'LW +0')
    f3_ax2.plot([flux_dict['mod_times'][tstp1], flux_dict['mod_times'][tstp1]], [flux_int_cnst_mol[tstp1], flux_int_cnst_mol[tstp1]], c='darkgray', linestyle='--')

    f3_ax2.scatter(flux_dict['mod_times'][tstp2], flux_int_cnst_mol[tstp2], s=150, c='r')
    f3_ax2.text(flux_dict['mod_times'][tstp2] + label_x_adj, flux_int_cnst_mol[tstp2], 'LW +3')
    f3_ax2.plot([flux_dict['mod_times'][tstp2], flux_dict['mod_times'][tstp2]], [flux_int_cnst_mol[tstp2], flux_int_cnst_mol[tstp2]], c='darkgray', linestyle='--')

    f3_ax2.scatter(flux_dict['mod_times'][tstp3], flux_int_cnst_mol[tstp3], s=150, c='r')
    f3_ax2.text(flux_dict['mod_times'][tstp3] + label_x_adj, flux_int_cnst_mol[tstp3], 'LW +6')
    f3_ax2.plot([flux_dict['mod_times'][tstp3], flux_dict['mod_times'][tstp3]], [flux_int_cnst_mol[tstp3], flux_int_cnst_mol[tstp3]], c='darkgray', linestyle='--')

    f3_ax2.scatter(flux_dict['mod_times'][tstp4], flux_int_cnst_mol[tstp4], s=150, c='r')
    f3_ax2.text(flux_dict['mod_times'][tstp4] + label_x_adj, flux_int_cnst_mol[tstp4], 'LW +9')
    f3_ax2.plot([flux_dict['mod_times'][tstp4], flux_dict['mod_times'][tstp4]], [flux_int_cnst_mol[tstp4], flux_int_cnst_mol[tstp4]], c='darkgray', linestyle='--')


    f3_ax2_txt = f3_ax1.text(-0.25, 1.0, '(b)', transform=f3_ax2.transAxes, fontsize='x-large')
    f3_ax2.set_ylabel('Regional CO$_{2}$ \n flux (mol h$^{-1}$)')
    f3_ax2.set_xlabel('Time (HH:MM)')
    f3_ax2.xaxis.set_major_formatter(HM_Fmt)

    
    f3_ax3 = fig3.add_subplot(gs[13:27, 0:22])
    f3_ax3.set_facecolor('olive')
    f3_ax3.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    f3_ax3.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(tot_flux_cnst[tstp1,:].shape), vmin=-1, vmax=1, cmap='bwr')
    f3_col =f3_ax3.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, flux_rate_cnst[tstp1,:], mask=~all_valid_c, zorder=2,  cmap=flux_int_cmap)
    f3_ax3.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=100, c='white', zorder=2)
    f3_ax3_txt_l4 = f3_ax3.annotate( 'L4',[tamar_fr_lon[l4_node_label_pos], tamar_fr_lat[l4_node_label_pos]], c='white')
    f3_ax3.set_title('LW+ 0 hrs',fontsize=30)
    f3_ax3_txt = f3_ax3.text(-0.13, 1.0, '(c)', transform=f3_ax3.transAxes, fontsize='x-large')


    area_label = (str(vaid_area) + ' km $^{2}$')
    f3_ax3_txt_area = f3_ax3.annotate( area_label,[tamar_fr_lon[text_node_label_pos], tamar_fr_lat[text_node_label_pos]],fontsize=20, c='white')


    f3_ax3.set_xlim([-4.28, -4.02])
    f3_ax3.set_ylim([50.23, 50.42])
    cax3 = f3_ax3.inset_axes([0.1, -0.1, 0.8, 0.03])
    f3_cb1 = fig3.colorbar(f3_col, cax=cax3, orientation='horizontal')
    f3_cb1.set_label(r'CO$_{2}$ flux (' 'mmol m$^{2}$ h$^{-1}$)')

    f3_ax3_xlab = f3_ax3.set_xlabel('$^{\circ}$ Lon')
    f3_ax3.xaxis.set_label_coords(-0.042, -0.02)
    f3_ax3_ylab = f3_ax3.set_ylabel('$^{\circ}$ Lat', rotation=0)
    f3_ax3.yaxis.set_label_coords(-0.08, 0.8)  

    f3_ax4 = fig3.add_subplot(gs[13:27, 28:50])
    f3_ax4.set_facecolor('olive')
    f3_ax4.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    f3_ax4.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(tot_flux_cnst[tstp2,:].shape), vmin=-1, vmax=1, cmap='bwr')
    f3_ax4.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, flux_rate_cnst[tstp2,:], mask=~all_valid_c, zorder=2, vmin=flux_min, vmax=flux_max, cmap=flux_int_cmap)
    f3_ax4.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=100, c='white', zorder=2)
    f3_ax4_txt_l4 = f3_ax4.annotate( 'L4',[tamar_fr_lon[l4_node_label_pos], tamar_fr_lat[l4_node_label_pos]], c='white')
    f3_ax4.set_title('LW+ 3 hrs',fontsize=30)
    f3_ax4_txt = f3_ax4.text(-0.13, 1.0, '(d)', transform=f3_ax4.transAxes, fontsize='x-large')
    area_label = (str(vaid_area) + ' km $^{2}$')
    f3_ax4_txt_area = f3_ax4.annotate( area_label,[tamar_fr_lon[text_node_label_pos], tamar_fr_lat[text_node_label_pos]],fontsize=20, c='white')

    f3_ax4.set_xlim([-4.28, -4.02])
    f3_ax4.set_ylim([50.23, 50.42])

    f3_ax4_xlab = f3_ax4.set_xlabel('$^{\circ}$ Lon')
    f3_ax4.xaxis.set_label_coords(-0.042, -0.02)
    f3_ax4_ylab = f3_ax4.set_ylabel('$^{\circ}$ Lat', rotation=0)
    f3_ax4.yaxis.set_label_coords(-0.08, 0.8)  
    cax4 = f3_ax4.inset_axes([0.1, -0.1, 0.8, 0.03])
    f3_cb2 = fig3.colorbar(f3_col, cax=cax4, orientation='horizontal')
    f3_cb2.set_label(r'CO$_{2}$ flux (' 'mmol m$^{2}$ h$^{-1}$)')

    f3_ax5 = fig3.add_subplot(gs[32:46, 0:22])
    f3_ax5.set_facecolor('olive')
    f3_ax5.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    flux_trp = f3_ax5.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(tot_flux_cnst[tstp3,:].shape), vmin=-1, vmax=1, cmap='bwr')
    f3_ax5.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, flux_rate_cnst[tstp3,:], mask=~all_valid_c, zorder=2, vmin=flux_min, vmax=flux_max, cmap=flux_int_cmap)
    f3_ax5.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=100, c='white', zorder=2)
    f3_ax5_txt_l4 = f3_ax5.annotate( 'L4',[tamar_fr_lon[l4_node_label_pos], tamar_fr_lat[l4_node_label_pos]], c='white')
    f3_ax5.set_title('LW+ 6 hrs',fontsize=30)
    f3_ax5_txt = f3_ax5.text(-0.13, 1.0, '(e)', transform=f3_ax5.transAxes, fontsize='x-large')
    area_label = (str(vaid_area) + ' km $^{2}$')
    f3_ax5_txt_area = f3_ax5.annotate( area_label,[tamar_fr_lon[text_node_label_pos], tamar_fr_lat[text_node_label_pos]],fontsize=20, c='white')

    f3_ax5.set_xlim([-4.28, -4.02])
    f3_ax5.set_ylim([50.23, 50.42])

    f3_ax5_xlab = f3_ax5.set_xlabel('$^{\circ}$ Lon')
    f3_ax5.xaxis.set_label_coords(-0.042, -0.02)
    f3_ax5_ylab = f3_ax5.set_ylabel('$^{\circ}$ Lat', rotation=0)
    f3_ax5.yaxis.set_label_coords(-0.08, 0.8)  
    cax5 = f3_ax5.inset_axes([0.1, -0.1, 0.8, 0.03])
    f3_cb3 = fig3.colorbar(f3_col, cax=cax5, orientation='horizontal')
    f3_cb3.set_label(r'CO$_{2}$ flux (' 'mmol m$^{2}$ h$^{-1}$)')

    f3_ax6 = fig3.add_subplot(gs[32:46, 28:50])
    f3_ax6.set_facecolor('olive')
    f3_ax6.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
    f3_ax6.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(tot_flux_cnst[tstp4,:].shape), vmin=-1, vmax=1, cmap='bwr')
    f3_col = f3_ax6.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, flux_rate_cnst[tstp4,:], mask=~all_valid_c, zorder=2, vmin=flux_min, vmax=flux_max, cmap=flux_int_cmap)
    f3_ax6.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=100, c='white', zorder=2)
    f3_ax6_txt_l4 = f3_ax6.annotate( 'L4',[tamar_fr_lon[l4_node_label_pos], tamar_fr_lat[l4_node_label_pos]], c='white')
    f3_ax6.set_title('LW+ 9 hrs',fontsize=30)
    f3_ax6_txt = f3_ax6.text(-0.13, 1.0, '(f)', transform=f3_ax6.transAxes, fontsize='x-large')
    area_label = (str(vaid_area) + ' km $^{2}$')
    f3_ax6_txt_area = f3_ax6.annotate( area_label,[tamar_fr_lon[text_node_label_pos], tamar_fr_lat[text_node_label_pos]],fontsize=20, c='white')

    f3_ax6.set_xlim([-4.28, -4.02])
    f3_ax6.set_ylim([50.23, 50.42])

    f3_ax6_xlab = f3_ax6.set_xlabel('$^{\circ}$ Lon')
    f3_ax6.xaxis.set_label_coords(-0.042, -0.02)
    f3_ax6_ylab = f3_ax6.set_ylabel('$^{\circ}$ Lat', rotation=0)
    f3_ax6.yaxis.set_label_coords(-0.08, 0.8)

    cax6 = f3_ax6.inset_axes([0.1, -0.1, 0.8, 0.03])
    f3_cb4 = fig3.colorbar(f3_col, cax=cax6, orientation='horizontal')
    f3_cb4.set_label(r'CO$_{2}$ flux (' 'mmol m$^{2}$ h$^{-1}$)')
    
    plt.tight_layout()

    fig3.savefig('tidal_flux_{}.png'.format(this_date.strftime('%Y-%m-%d')), dpi=180)
    if save_pdf:
        fig3.savefig('tidal_flux_{}.pdf'.format(this_date.strftime('%Y-%m-%d')), dpi=60)
    plt.close() 

fluxsave[3,2]/fluxsave[3,0]







### Landshcutzer and long term flux plots


#this is with fixed wind
const4 = np.load('new_out_data_all_7_1.npy', allow_pickle=True).item()

time_dt4 = const4['time_dt']
flux_int4 = const4['flux_int']/1e6 #this is umol hr-1 integrated divide by 1e6 to get moles
valid_const4 = const4['all_valid']
fco2_const4 = const4['mean_fc02']
stdfco2_const4 = const4['std_fc02']

all_valid_nodes_const4 = np.where(valid_const4)[0]
all_valid_ele_const4 = np.where(np.isin(tamar_fr_triangles,all_valid_nodes_const4))[0]

all_valid_c_const4 = np.zeros(len(tamar_fr_lonc), dtype=bool)
all_valid_c_const4[all_valid_ele_const4] = True



#loop through flux and get weekly variables
hrs_week=24*7
start=0
loop_it=np.floor(len(time_dt4)/hrs_week)
time_dt4_week=[0] * int(loop_it)
flux_int4_week=[0] * int(loop_it)
for v in range(0,int(loop_it)):
    hours_select = list(range((hrs_week*(v-1))+1,(hrs_week*v)))
    hour_int=np.array(hours_select)
    flux_int4_week[v-1]=np.average(flux_int4[hour_int], axis=0)
    week_dt_stamp=[hours_select[1]+(hrs_week/2)]
    y=int(week_dt_stamp[0])
    time_dt4_week[v-1]=time_dt4[y]
flux_int4_week_floar=np.array(flux_int4_week)



#this is with varying wind
var = np.load('new_out_data_all_var.npy', allow_pickle=True).item()

time_dt_var = var['time_dt']
flux_int_var = var['flux_int']/1e6#this is umol hr-1 integrated divide by 1e6 to get moles

cut_off = 10**18
remove = np.logical_or(flux_int_var > cut_off, flux_int_var < -cut_off)

cut_off = 10**18
remove2 = np.logical_or(flux_int_var > cut_off, flux_int_var < -cut_off)


# this is landshutzer! with non varying winds

land_ls = np.load('out_data_landsch_7_1_fmolhr.npy', allow_pickle=True).item()
time_land = land_ls['time_dt']
flux_land = land_ls['flux_int']/1e6#this is umol hr-1 integrated divide by 1e6 to get moles
valid_land = land_ls['all_valid']


all_valid_nodes_land = np.where(valid_land)[0]
all_valid_ele_land = np.where(np.isin(tamar_fr_triangles,all_valid_nodes_land))[0]

all_valid_c_land = np.zeros(len(tamar_fr_lonc), dtype=bool)
all_valid_c_land[all_valid_ele_land] = True


#loop through flux and get weekly variables
hrs_week=24*7
start=0
loop_it=np.floor(len(time_land)/hrs_week)
time_land_week=[0] * int(loop_it)
flux_land_week=[0] * int(loop_it)
for v in range(0,int(loop_it)):
    hours_select = list(range((hrs_week*(v-1))+1,(hrs_week*v)))
    hour_int_land=np.array(hours_select)
    flux_land_week[v-1]=np.average(flux_land[hour_int_land], axis=0)
    week_dt_stamp=[hours_select[1]+(hrs_week/2)]
    y=int(week_dt_stamp[0])
    time_land_week[v-1]=time_land[y]
flux_land_week_floar=np.array(flux_land_week)




#this is the flux calculated with landshcutzer and l4 measurements
#scaled up by area
land_l4_data = np.load('pt_flux_fmolh.npy', allow_pickle=True).item()
scaled_flux_dt = land_l4_data['time_dt']
l4_flux_scaled = land_l4_data['obs_cst_flux']/1e6#this is umol hr-1 integrated divide by 1e6 to get moles
land_flux_scaled = land_l4_data['landsch_cst_flux']/1e6#this is umol hr-1 integrated divide by 1e6 to get moles
hrs_week=24*7
start=0
loop_it=np.floor(len(l4_flux_scaled)/hrs_week)
scaled_flux_dt_week=[0] * int(loop_it)
l4_flux_scaled_week=[0] * int(loop_it)
land_flux_scaled_week=[0] * int(loop_it)

for v in range(0,int(loop_it)):
    hours_select = list(range((hrs_week*(v-1))+1,(hrs_week*v)))
    
    hour_int_scaled_flux=np.array(hours_select)
    l4_flux_scaled_week[v-1]=np.average(l4_flux_scaled[hour_int_scaled_flux], axis=0)
    
    hour_int_scaled_flux=np.array(hours_select)
    land_flux_scaled_week[v-1]=np.average(land_flux_scaled[hour_int_scaled_flux], axis=0)
    
    week_dt_stamp_scaled=[hours_select[1]+(hrs_week/2)]
    y=int(week_dt_stamp_scaled[0])
    scaled_flux_dt_week[v-1]=scaled_flux_dt[y]
    
l4_flux_scaled_week_floar=np.array(l4_flux_scaled_week)
l4_flux_scaled_week_floar_mol=np.array(l4_flux_scaled_week)

land_flux_scaled_week_floar=np.array(land_flux_scaled_week)
land_flux_scaled_week_floar_mol=np.array(land_flux_scaled_week)




wrf_l4_data = np.load('wrf_l4_wnd.npy', allow_pickle=True)
buoy_dt = wrf_l4_data[0]
buoy_dt, buoy_ind = np.unique(buoy_dt, return_index=True)

wind_spd = np.asarray(wrf_l4_data[1], dtype=float)[buoy_ind]

l4_fC02_file = 'Quest_2016_fCO2_L4_only.csv'


fC02_data = np.loadtxt(l4_fC02_file, delimiter=',', skiprows=1)
fC02_dt_bst = np.asarray([dt.datetime(int(this_row[0]), int(this_row[1]), int(this_row[2]), int(this_row[3]), int(this_row[4]), int(this_row[5])) for this_row in fC02_data])
fC02_dt_raw = np.asarray([this_dt - dt.timedelta(hours=1) for this_dt in fC02_dt_bst])
fC02_day = np.asarray([dt.datetime(2016, this_dt.month, this_dt.day) for this_dt in fC02_dt_raw])
fC02_day_u = np.unique(fC02_day)

fC02_dt = []
fC02 = []
fC02_atmos = []
pC02 = []
pC02_atmos = []
fC02_lon = []
fC02_lat = []

for this_day in fC02_day_u:
    choose_dt = fC02_day == this_day

    poss_dates = fC02_dt_raw[choose_dt]
    fC02_dt.append(poss_dates[int(np.floor(len(poss_dates)/2))])

    fC02_atmos.append(np.mean(fC02_data[choose_dt,9]))
    fC02.append(np.mean(fC02_data[choose_dt,8]))
    pC02.append(np.mean(fC02_data[choose_dt,6]))
    pC02_atmos.append(np.mean(fC02_data[choose_dt,7]))

    fC02_lon.append(np.mean(fC02_data[choose_dt,13]))
    fC02_lat.append(np.mean(fC02_data[choose_dt,12]))

fC02_dt = np.asarray(fC02_dt)
fC02 = np.asarray(fC02)
fC02_atmos = np.asarray(fC02_atmos)
pC02 = np.asarray(pC02)
pC02_atmos = np.asarray(pC02_atmos)
fC02_lon = np.asarray(fC02_lon)
fC02_lat = np.asarray(fC02_lat)

ls_file = 'MPI-ULB-SOM_FFN_clim.nc'
ls_nc = nc.Dataset(ls_file)

bounds = [[-4.3, -4], [50.1, 50.2]]

ls_lon = ls_nc.variables['lon'][:]
ls_lat = ls_nc.variables['lat'][:]

choose_lon = np.logical_and(ls_lon >= bounds[0][0], ls_lon <= bounds[0][1])
choose_lat = np.logical_and(ls_lat >= bounds[1][0], ls_lat <= bounds[1][1])

lon_c, lat_c = np.meshgrid(ls_lon[choose_lon], ls_lat[choose_lat])

landshc_pc02 = np.squeeze(ls_nc['pco2'][:,choose_lat, choose_lon])
landshc_bound_lon = ls_nc['lon_bnds'][choose_lon]
landshc_bound_lat = ls_nc['lat_bnds'][choose_lat]

monthdays = np.asarray([dt.datetime(2016,m,15) for m in np.arange(1,13)])


plt.rcParams.update({'font.size': 26})
Month_Fmt = mdates.DateFormatter('%B')

# fig3 = plt.figure(figsize=(24,30))
# gs = fig3.add_gridspec(3, 1)

# fig3.subplots_adjust(top=0.98)
# fig3.subplots_adjust(left=0.1)
# fig3.subplots_adjust(right=0.9)
# fig3.subplots_adjust(wspace=0.15)
# fig3.subplots_adjust(hspace=0.15)
# fig3.subplots_adjust(bottom=0.08)

# f3_ax1 = fig3.add_subplot(gs[0, :])
# f3_ax1.plot(monthdays, landshc_pc02, c='C0', linewidth=2)
# f3_ax1.plot(fC02_dt, pC02, c='C1', linewidth=2)
# f3_ax1.plot(fC02_dt, pC02_atmos, c='C2', linewidth=2)
# f3_ax1.set_xlim([np.min(time_dt4), np.max(time_dt4)])
# f3_ax1.xaxis.set_major_formatter(Month_Fmt)
# f3_ax1.set_ylabel('pC02, $\mu atm$')
# f3_ax1_txt = f3_ax1.text(0.03, 0.9, '(A)', transform=f3_ax1.transAxes, fontsize='x-large')


# f3_ax2 = fig3.add_subplot(gs[1, :])
# #f3_ax3.plot(buoy_dt[~remove_ex], wind_spd[~remove_ex])

# f3_ax2_2 = f3_ax2.twinx()
# f3_ax2.plot(fC02_dt, pC02, c='C1', linewidth=2)  # Changed from fC02 to pC02 to be same as atmospheric
# #f3_ax2.scatter(fC02_dt, pC02, s=150, c='C1')
# f3_ax2.plot(fC02_dt, pC02_atmos, c='C2', linewidth=2)
# f3_ax2_2.plot(buoy_dt, wind_spd,c='C4', linewidth=2)
# f3_ax2.xaxis.set_major_formatter(Month_Fmt)
# f3_ax2.set_xlim([np.min(time_dt4), np.max(time_dt4)])
# f3_ax2.set_ylabel('pC02, $\mu atm$')
# f3_ax2_2.set_ylabel('Wind speed, ms$^{-1}$')
# f3_ax2_txt = f3_ax2.text(0.03, 0.9, '(B)', transform=f3_ax2.transAxes, fontsize='x-large')

# f3_ax3 = fig3.add_subplot(gs[2:, :])
# #f3_ax1.plot(time_dt, flux_int)
# f3_ax3.plot(time_dt4, np.zeros(len(time_dt4)), c='k', linestyle='--')
# f3_ax3.plot(time_dt4, -flux_int4, c='darkgreen')
# f3_ax3.plot(time_dt_var[~remove2], -flux_int_var[~remove2], c='r')
# f3_ax3.set_xlim([np.min(time_dt4), np.max(time_dt4)])
# f3_ax3.xaxis.set_major_formatter(Month_Fmt)
# f3_ax3.set_ylabel('Integrated C02 flux, mmol h$^{-1}$')
# f3_ax3_txt = f3_ax3.text(0.03, 0.2, '(C)', transform=f3_ax3.transAxes, fontsize='x-large')

# fig3.savefig('longterm_plot.png', dpi=180)
# if save_pdf:
#     fig3.savefig('longterm_plot.pdf', dpi=60)
# plt.close() 

#valid long term in m2
vaid_area_longterm =(np.sum(tamar_fr_art1[all_valid_nodes_const4], axis=0))#convert from m2 to km 2



fig6 = plt.figure(figsize=(24,15))
gs = fig6.add_gridspec(1, 1)

f6_ax1 = fig6.add_subplot(gs[0, :])
f6_ax1_2 = f6_ax1.twinx()
plt.axvline(dt.datetime(2016,6, 10),color='k', linestyle='dashed')
plt.axvline(dt.datetime(2016,9, 21),color='k', linestyle='dashed')
lns1=f6_ax1.plot(time_dt4, fco2_const4,'-o',  color='C0', linewidth=2, markersize=4, label = 'This study')
#this plots the range- removed
#error_fill=f6_ax1.fill_between(time_dt4, fco2_const4-(stdfco2_const4*2), fco2_const4+(stdfco2_const4*2), alpha=0.2, color='C0')

error_fill=f6_ax1.fill_between(time_dt4, fco2_const4-(21.88 ), fco2_const4+(21.88 ), alpha=0.2, color='C0')
#note added time range to pco2 plots to stop them plotting off the plot when joined with lines
lns2=f6_ax1.plot(fC02_dt[1:-3], fC02[1:-3],'-o', c='green', linewidth=2, markersize=12,label = 'L4')
lns3=f6_ax1_2.plot(monthdays[2:10], landshc_pc02[2:10],'-o', c='orange', linewidth=2, markersize=12,label ='Landschützer coastal product' )
lns4=f6_ax1.plot(fC02_dt[1:-3], fC02_atmos[1:-3],'-s', c='grey', linewidth=2, markersize=12,label = 'L4 atmosphere')
f6_ax1.set_xlim([np.min(time_dt4), np.max(time_dt4)])
f6_ax1.xaxis.set_major_formatter(Month_Fmt)
f6_ax1.set_ylabel(r'$fCO_{2}$' ' (' '\u03BC' 'atm)',fontsize=30)
#f6_ax1.text(-0.1, 1.01, '(a)',transform=f6_ax1.transAxes,fontsize=50)
f6_ax1.set_xlabel('Time',fontsize=30)



# added these three lines
lns = lns1+lns2+lns3+lns4
labs = [l.get_label() for l in lns]
f6_ax1.legend(lns, labs, loc=[0.45,0.76])
f6_ax1_2.set_ylabel(r'$pCO_{2}$' ' (' '\u03BC' 'atm)',fontsize=30)
f6_ax1.set_ylim([330, 525])
f6_ax1_2.set_ylim([330, 525])



# #comment start on 12b 

# f6_ax2 = fig6.add_subplot(gs[1:, :])
# #f3_ax1.plot(time_dt, flux_int)
# # f6_ax2.plot(time_dt4, -flux_int4, c='darkgreen')

# #this is hourly model run averaged to weekly with monthly winds
# f6_ax2.plot(time_dt4_week, -flux_int4_week_floar,'-o',  c='C0',linewidth=2, markersize=18)

# # #this is hourly model run averaged to weekly with monthly winds - need to update to replot!
# # f6_ax2.plot(time_land_week, -flux_land_week_floar,'-o',  c='red', markersize=18)

# #not spatially varying but using land and L4 CO2 scaled
# #l4
# f6_ax2.plot(scaled_flux_dt_week, -l4_flux_scaled_week_floar_mol,'-o', linewidth=2, c='green', markersize=18)
# #landshutzer
# f6_ax2.plot(scaled_flux_dt_week, -land_flux_scaled_week_floar_mol,'-o',linewidth=2,  c='orange', markersize=18)
# f6_ax2.plot(time_dt4, np.zeros(len(time_dt4)), c='k', linestyle='--')
# f6_ax2.set_xlim([np.min(time_dt4), np.max(time_dt4)])
# f6_ax2.xaxis.set_major_formatter(Month_Fmt)
# f6_ax2.set_ylabel('Regional CO$_{2}$ \n flux (''mol h$^{-1}$)',fontsize=30)
# f6_ax2.set_xlabel('Time',fontsize=30)
# f6_ax2.text(-0.1, 1.01, '(b)',transform=f6_ax2.transAxes,fontsize=50)
# f6_ax2.legend(['This study','L4','Landschützer coastal product'])
# # f6_ax2.set_ylim([-50000, 10000])

# # f6_ax2_2 = f6_ax2.twinx()

# #comment end 12b

fig6.savefig('longterm_flux_compared.png', dpi=180)
if save_pdf:
    fig6.savefig('longterm_flux_compared.pdf', dpi=60)
plt.close() 




# this plot is just to check valid area

fig3 = plt.figure(figsize=(22,20), constrained_layout=True)
gs = fig3.add_gridspec(1, 2)
f3_ax2 = fig3.add_subplot(gs[0, 1])
f3_ax2.set_facecolor('olive')
f3_ax2.triplot(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, c='lightgray')
f3_ax2.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, np.zeros(sal.shape), vmin=-1, vmax=1, cmap='bwr')
f3_col2 = f3_ax2.tripcolor(tamar_fr_lon, tamar_fr_lat, tamar_fr_triangles, fc02, mask=~all_valid_c_const4, zorder=2, cmap=flux_int_cmap)
f3_ax2.scatter(tamar_fr_lon[l4_node], tamar_fr_lat[l4_node], s=150, c='white', zorder=2)
f3_ax2_txt_l4 = f3_ax2.annotate( 'L4',[tamar_fr_lon[l4_node_label_pos], tamar_fr_lat[l4_node_label_pos]], c='white')
f3_ax2_txt = f3_ax2.text(0.00, 1.03, '(b)', transform=f3_ax2.transAxes, fontsize='x-large')
f3_ax2_txt.set_bbox(dict(facecolor='white', edgecolor='white'))
f3_ax2.set_xlim([-4.28, -4.02])
f3_ax2.set_ylim([50.23, 50.4])
f3_ax2_xlab = f3_ax2.set_xlabel('$^{\circ}$ Lon')
f3_ax2.xaxis.set_label_coords(-0.05, -0.02)
f3_ax2_ylab = f3_ax2.set_ylabel('$^{\circ}$ Lat', rotation=0)
f3_ax2.yaxis.set_label_coords(-0.08, 0.9)  

cax2 = f3_ax2.inset_axes([0.1, -0.1, 0.8, 0.03])
f3_cb2 = fig3.colorbar(f3_col2, cax=cax2, orientation='horizontal')
f3_cb2.set_label(r'$fCO_{2 (sw)}$' '(' '\u03BC' 'atm)')
fig3.savefig('longterm_flux_valid_area.png', dpi=180)

vaid_area_longterm = round((np.sum(tamar_fr_art1[all_valid_nodes_const4], axis=0))/1e6)#convert from m2 to km 2
