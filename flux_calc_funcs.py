import numpy as np

import PyFVCOM as pf



def flux_c02_calc(u10, T, deltaC02, S):
    Tk = T + 273.15

    k660=(0.222*(u10**2) +0.333*u10) # from nightingale 2000
    Schcw=(2073.1 -(125.62*T) + (3.6276*(T**2)) - (0.04321*(T**3))) #from wannikhoif 1992 - for seawater
    Schdep=(Schcw/660)**(-0.5) #johnson 2010 - Dimensionless
    Kw=Schdep*k660 #k with units of cm hr-1
    k0=np.exp(-60.2409 + 93.4517*(100/Tk) + 23.3585*np.log(Tk/100) +S*(0.023517 -0.023656*(Tk/100) + 0.0047036*((Tk/100)**2))) #K0 with units(mol l-1 atm-1) weiss 1974
    scal=(12*1000)/(10**6) # scaling K0- sclaing factor to convert from mol L-1 atm-1 to gC m_3 micro atm-1
    scal2=(1/100) #scaling Kw from cm hr-1 to m hr-1
 
    #transfer coefficent
    TR=Kw*k0*scal*scal2 #Units of gC m-2 hr-1 atm-1
    TRmo=TR*24*365/12   #Units of gC m-2 month-1 atm-1

    #calculate flux using delta co2
    Fhr=TR*deltaC02 #Units of gC m-2 hr-1
    Fmo=TRmo*deltaC02 #Units of gC m-2 mo-1
 
    # Convert the flux of CO2 to umol m-2 hr-1
    Fmolhr=Fhr*(1/12)*(10**6)

    return Fhr, Fmo, Fmolhr


def get_mod_flux(mod_fr, calc_time, deltaC02_node, wind_speed=None):
    
    time_ind = mod_fr.closest_time(calc_time)
   
    T = mod_fr.data.temp[time_ind,0,:]
    S = mod_fr.data.salinity[time_ind,0,:]

    if wind_speed == None:
        wind_spd_ele = np.sqrt(mod_fr.data.uwind_speed[:,0,:]**2 + mod_fr.data.vwind_speed[:,0,:]**2)
        
        close_eles = wind_fr.closest_element([wind_fr.grid.lon, wind_fr.grid.lat])        
        wind_spd_node = wind_spd_ele[:,close_eles]
    else:
        wind_spd_node = wind_speed * np.ones(T.shape)

    grid_flux = flux_c02_calc(wind_spd_node, T, deltaC02_node, S)

    return grid_flux




