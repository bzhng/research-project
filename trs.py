import os
from time import time
import multiprocessing
import netCDF4 as nc
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt



rho = 1025
g = 9.81
T = (12*60+25)*60
omega = np.pi*2/T 
phi = np.pi/4
A = 1
Cp = 1

epoch = np.datetime64('1970-01-01T00:00:00')



    
    
   
    
    
    
    
def peak_trough_idx(height,t,latitude,longitude):
    '''
    indices of the peaks and troughs
    
    '''
#     count = 0
    lat_error = []
    lon_error = []
    
    
    secs = (t-epoch)/1e9 # dtype is timedelta[s] not timedelta[ns]
    secs = secs.astype('float64')
    
    step = 50
    istep = 35
    N = (secs[-1]-secs[0])/T
    
    if N%1 > 0.4 and N%1 < 0.9:
        
        n = int(N)
            
        idx_troughs = np.zeros((n+1,len(latitude),len(longitude)))
        idx_peaks = np.zeros((n+1,len(latitude),len(longitude)))
        idx_troughs[:] = np.nan
        idx_peaks[:] = np.nan

        
    else:
    
        idx_troughs = np.zeros((round(N),len(latitude),len(longitude)))
        idx_peaks = np.zeros((round(N),len(latitude),len(longitude)))
        idx_troughs[:] = np.nan
        idx_peaks[:] = np.nan

    
    
    
    for j, lat in enumerate(latitude):
        for k, long in enumerate(longitude):
            
            if np.isnan(height[0,j,k]) == True:
                pass
            
            
            
            if height[0,j,k]>height[1,j,k] or height[0,j,k]==height[1,j,k] and height[1,j,k]>height[2,j,k]:
                
                if N%1 > 0.4 and N%1 < 0.9:
                    peak_shape = n
                    trough_shape = n+1
                    
                    trough = min(height[:istep,j,k])
                    idx_trough = [x for x, z in enumerate(height[:istep,j,k]) if z == trough][0]
                    if height[idx_trough,j,k] == height[istep-1,j,k]:
                        trough = min(height[:istep+step,j,k])
                        idx_trough = [x for x, z in enumerate(height[:istep+step,j,k]) if z == trough][0]
                    idx_troughs[0,j,k] = idx_trough

                    
                
                    for i in range(1,trough_shape):
                        
                        peak_slice = height[idx_trough:idx_trough+step,j,k]
                        peak = max(peak_slice)
                        idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        if height[idx_peak,j,k] == peak_slice[-1]:
                            peak_slice = height[idx_trough:idx_trough+step+step,j,k]
                            peak = max(peak_slice)
                            idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        idx_peaks[i-1,j,k] = idx_peak

                       
                        trough_slice = height[idx_peak:idx_peak+step,j,k]
                        trough = min(trough_slice)
                        idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        if height[idx_trough,j,k] == trough_slice[-1]:
                            trough_slice = height[idx_peak:idx_peak+step+step,j,k]
                            trough = min(trough_slice)
                            idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        idx_troughs[i,j,k] = idx_trough

                        
                       
                    
                
                else:
                    peak_shape = round(N)
                    trough_shape = round(N)
                
                    trough = min(height[:istep,j,k])
                    idx_trough = [x for x, z in enumerate(height[:istep,j,k]) if z == trough][0]
                    if height[idx_trough,j,k] == height[istep-1,j,k]:
                        trough = min(height[:istep+step,j,k])
                        idx_trough = [x for x, z in enumerate(height[:istep+step,j,k]) if z == trough][0]
                    idx_troughs[0,j,k] = idx_trough

                    
               
                
                    for i in range(1,trough_shape):

                        
                        peak_slice = height[idx_trough:idx_trough+step,j,k]
                        peak = max(peak_slice)
                        idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        if height[idx_peak,j,k] == peak_slice[-1]:
                            peak_slice = height[idx_trough:idx_trough+step+step,j,k]
                            peak = max(peak_slice)
                            idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        idx_peaks[i-1,j,k] = idx_peak
                        
                        trough_slice = height[idx_peak:idx_peak+step,j,k]
                        trough = min(trough_slice)
                        idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        if height[idx_trough,j,k] == trough_slice[-1]:
                            trough_slice = height[idx_peak:idx_peak+step+step,j,k]
                            trough = min(trough_slice)
                            idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        idx_troughs[i,j,k] = idx_trough

                       
                        
                     
                    peak_slice = height[idx_trough:idx_trough+step,j,k]
                    peak = max(peak_slice)
                    idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                    if height[idx_peak,j,k] == peak_slice[-1]:
                        peak_slice = height[idx_trough:idx_trough+step+step,j,k]
                        peak = max(peak_slice)
                        idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                    idx_peaks[-1,j,k] = idx_peak
                    
                    
                
                
                
            elif height[0,j,k]<height[1,j,k] or height[0,j,k]==height[1,j,k] and height[1,j,k]<height[2,j,k]:
                
                if N%1 > 0.4 and N%1 < 0.9:
                    peak_shape = n+1
                    trough_shape = n
                    
                    peak = max(height[:istep,j,k])
                    idx_peak = [x for x, z in enumerate(height[:istep,j,k]) if z == peak][0]
                    if height[idx_peak,j,k] == height[istep-1,j,k]:
                        peak = max(height[:istep+step,j,k])
                        idx_peak = [x for x, z in enumerate(height[:istep+step,j,k]) if z == peak][0]
                    idx_peaks[0,j,k] = idx_peak

                    


                    for i in range(1,peak_shape):

                        trough_slice = height[idx_peak:idx_peak+step,j,k]
                        trough = min(trough_slice)
                        idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        if height[idx_trough,j,k] == trough_slice[-1]:
                            trough_slice = height[idx_peak:idx_peak+step+step,j,k]
                            trough = min(trough_slice)
                            idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        idx_troughs[i-1,j,k] = idx_trough
                        
                        
                        peak_slice = height[idx_trough:idx_trough+step,j,k]
                        peak = max(peak_slice)
                        idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        if height[idx_peak,j,k] == peak_slice[-1]:
                            peak_slice = height[idx_trough:idx_trough+step+step,j,k]
                            peak = max(peak_slice)
                            idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        idx_peaks[i,j,k] = idx_peak
                        
                        

                 
                    
                
                else:
                    peak_shape = round(N)
                    trough_shape = round(N)
                
                    peak = max(height[:istep,j,k])
                    idx_peak = [x for x, z in enumerate(height[:istep,j,k]) if z == peak][0]
                    if height[idx_peak,j,k] == height[istep-1,j,k]:
                        peak = max(height[:istep+step,j,k])
                        idx_peak = [x for x, z in enumerate(height[:istep+step,j,k]) if z == peak][0]
                    idx_peaks[0,j,k] = idx_peak


                    for i in range(1,peak_shape):

                        
                        trough_slice = height[idx_peak:idx_peak+step,j,k]
                        trough = min(trough_slice)
                        idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        if height[idx_trough,j,k] == trough_slice[-1]:
                            trough_slice = height[idx_peak:idx_peak+step+step,j,k]
                            trough = min(trough_slice)
                            idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                        idx_troughs[i-1,j,k] = idx_trough
                        
                        
                        peak_slice = height[idx_trough:idx_trough+step,j,k]
                        peak = max(peak_slice)
                        idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        if height[idx_peak,j,k] == peak_slice[-1]:
                            peak_slice = height[idx_trough:idx_trough+step+step,j,k]
                            peak = max(peak_slice)
                            idx_peak = [x for x, z in enumerate(peak_slice) if z == peak][0] + idx_trough
                        idx_peaks[i,j,k] = idx_peak
                        
                        
                     
                    trough_slice = height[idx_peak:idx_peak+step,j,k]
                    trough = min(trough_slice)
                    idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                    if height[idx_trough,j,k] == trough_slice[-1]:
                        trough_slice = height[idx_peak:idx_peak+step+step,j,k]
                        trough = min(trough_slice)
                        idx_trough = [x for x, z in enumerate(trough_slice) if z == trough][0] + idx_peak
                    idx_troughs[-1,j,k] = idx_trough
                    
                        

                        
      
            peak_filter = idx_peaks[idx_peaks[:,j,k]<len(height[:,j,k]),j,k]
            trough_filter = idx_troughs[idx_troughs[:,j,k]<len(height[:,j,k]),j,k]
            
            if len(peak_filter)==0 and len(trough_filter)==0:
                pass
            
            elif len(peak_filter)!=len(np.unique(peak_filter)) or len(trough_filter)!=len(np.unique(trough_filter)) or peak_filter[-1]==trough_filter[-1]:
#                 count += 1
                lat_error.append(j)
                lon_error.append(k)
#                 print(f'{count} ')
#                 print(f'lat={lat}, index={j}; \nlon={long}, index={k}\n\n')
                   
                idx_peaks[len(peak_filter[peak_filter<peak_filter[-1]]):,j,k] = np.nan
                idx_troughs[len(trough_filter[trough_filter<trough_filter[-1]]):,j,k] = np.nan
        
            
                
            
    strange = np.array([lat_error,lon_error])
                
                   
                
               

                    
    return idx_peaks,idx_troughs,strange
   
    
    
    
    
def write_peak_trough(idx_peaks,idx_troughs,latitude,longitude,filename):
    '''
    
    '''    
        
    indices = xr.Dataset({
        'peak': xr.DataArray(
            data=idx_peaks,
            coords={'p_occur': np.arange(1,idx_peaks.shape[0]+1),
                    'lat': latitude,
                    'lon': longitude},
            dims=['p_occur','lat','lon']),
        
        'trough': xr.DataArray(
            data=idx_troughs,
            coords={'t_occur': np.arange(1,idx_troughs.shape[0]+1),
                    'lat': latitude,
                    'lon': longitude},
            dims=['t_occur','lat','lon'])})
                                   

    indices.to_netcdf(path=f'data/processed/{filename}')
    
    return










def R(height,idx_peaks,idx_troughs,peak_occur,trough_occur,latitude,longitude):
    '''
    
    '''
    occur = len(peak_occur) # shape of peak_occur=trough_occur
    

    R_map = np.zeros((occur,2,len(latitude),len(longitude))) 
    R_map[:] = np.nan

    for j, lat in enumerate(latitude):
        for k, long in enumerate(longitude):
            
            if np.isnan(height[0,j,k]) == True:
                pass
            if np.isnan(idx_peaks[0,j,k]) == True:
                pass
            if np.isnan(idx_troughs[0,j,k]) == True:
                pass
            
            
            
            else:
                
                peaks = idx_peaks[idx_peaks[:,j,k]<=len(height[:,j,k]),j,k].astype(int)
                troughs = idx_troughs[idx_troughs[:,j,k]<=len(height[:,j,k]),j,k].astype(int)

                R_array = np.zeros((occur,2)) 
                R_array[:] = np.nan

                if len(peaks) == len(troughs):

                    if peaks[0] > troughs[0]:
                        for i in range(len(peaks)-1):
                            R_array[i] = height[peaks[i],j,k]-height[troughs[i],j,k],height[peaks[i],j,k]-height[troughs[i+1],j,k]
                        R_array[-1,0] = height[peaks[-1],j,k] - height[troughs[-1],j,k]


                    if peaks[0] < troughs[0]:
                        R_array[0,1] = height[peaks[0],j,k] - height[troughs[0],j,k]
                        for i in range(1,len(peaks)):
                            R_array[i] = height[peaks[i],j,k]-height[troughs[i-1],j,k],height[peaks[i],j,k]-height[troughs[i],j,k]        



                if len(peaks) == len(troughs)+1:
                    R_array[0,1] = height[peaks[0],j,k] - height[troughs[0],j,k]
                    for i in range(1,len(peaks)-1):
                            R_array[i] = height[peaks[i],j,k]-height[troughs[i-1],j,k],height[peaks[i],j,k]-height[troughs[i],j,k]        
                    R_array[-1,0] = height[peaks[-1],j,k] - height[troughs[-1],j,k]


                if len(peaks)+1 == len(troughs):
                    for i in range(len(peaks)):
                        R_array[i] = height[peaks[i],j,k]-height[troughs[i],j,k],height[peaks[i],j,k]-height[troughs[i+1],j,k]


                R_map[:,:,j,k] = R_array
                
        
    return R_map









def P_in(R,t,gamma):
    '''
    instantaneous power density
    '''
    return rho*g*(R/2)*(np.cos(phi)*np.cos(omega*t-phi-gamma)-np.cos(omega*t-gamma))*(R/2)*omega*np.cos(phi)*np.sin(omega*t-phi-gamma)

    ### phase dependence of actual elevation data





def PD_t(range_array,dt,gamma,idx_peaks,idx_troughs,peak_occur,trough_occur,latitude,longitude): 
    '''
    '''
    
    secs = (dt-epoch)/1e9 # dtype is timedelta[s] not timedelta[ns]
    secs = secs.astype('float64')
    secs = secs - secs[0] ### shift time to start at zero
    
    
    # requires removal/filtering of trailing nan if the first point at each lat,long is not nan
    pd = np.zeros((len(secs),len(latitude),len(longitude)))
    pd[:] = np.nan
    
    # initialised with array of epoch datetimes, requires removal of trailing epochs
    t_pd = np.zeros((len(secs),len(latitude),len(longitude)),dtype='datetime64[s]') 
    
    


    for j, lat in enumerate(latitude):
        for k, long in enumerate(longitude):

#             print(j,k)
            
            if np.isnan(idx_peaks[0,j,k]) == True:
                pass
            if np.isnan(idx_troughs[0,j,k]) == True:
                pass

            else:
                occur = len(peak_occur) # shape of the two 'occur(ence)' arrays are the same
                peaks = idx_peaks[idx_peaks[:,j,k]<=len(secs),j,k].astype(int)
                troughs = idx_troughs[idx_troughs[:,j,k]<=len(secs),j,k].astype(int)
                
                

                power_density = []

                if len(peaks) == len(troughs):

                    if peaks[0] > troughs[0]:
                        if np.isnan(range_array[0,0,j,k]) == True: 
                            pass
                        else:  
                            for i in range(len(peaks)-1):
                                t = secs[troughs[i]:peaks[i]]
                                power = P_in(range_array[i,0,j,k],t,gamma[j,k])
                                for p in power:
                                    power_density.append(p)
                                t = secs[peaks[i]:troughs[i+1]]
                                power = P_in(range_array[i,1,j,k],t,gamma[j,k])
                                for p in power:
                                    power_density.append(p)

                            t = secs[troughs[-1]:peaks[-1]+1]
                            power = P_in(range_array[i+1,0,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)
                                
                            t_dt = dt[troughs[0]:peaks[-1]+1]

                    if peaks[0] < troughs[0]:
                        if np.isnan(range_array[0,1,j,k]) == True:
                            pass
                        else:
                            t = secs[peaks[0]:troughs[0]]
                            power = P_in(range_array[0,1,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)
                            for i in range(1,len(peaks)-1):
                                t = secs[troughs[i-1]:peaks[i]]
                                power = P_in(range_array[i,0,j,k],t,gamma[j,k])
                                for p in power:
                                    power_density.append(p)
                                t = secs[peaks[i]:troughs[i]]
                                power = P_in(range_array[i,1,j,k],t,gamma[j,k])
                                for p in power:
                                    power_density.append(p)

                            t = secs[troughs[-2]:peaks[-1]]
                            power = P_in(range_array[i+1,0,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)
                            t = secs[peaks[-1]:troughs[-1]+1]
                            power = P_in(range_array[i+1,1,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)
                                
                            t_dt = dt[peaks[0]:troughs[-1]+1]


                if len(peaks) == len(troughs)+1:
                    if np.isnan(range_array[0,1,j,k]) == True:
                        pass
                    else:
                        t = secs[peaks[0]:troughs[0]]
                        power = P_in(range_array[0,1,j,k],t,gamma[j,k])
                        for p in power:
                            power_density.append(p)
                        for i in range(1,len(peaks)-1):
                            t = secs[troughs[i-1]:peaks[i]]
                            power = P_in(range_array[i,0,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)
                            t = secs[peaks[i]:troughs[i]]
                            power = P_in(range_array[i,1,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)
                        t = secs[troughs[-1]:peaks[-1]+1]
                        power = P_in(range_array[i+1,0,j,k],t,gamma[j,k])
                        for p in power:
                            power_density.append(p)

                        t_dt = dt[peaks[0]:peaks[-1]+1]




                if len(peaks)+1 == len(troughs):
                    if np.isnan(range_array[0,0,j,k]) == True:
                        pass
                    else:
                        for i in range(len(peaks)-1):
                            t = secs[troughs[i]:peaks[i]]
                            power = P_in(range_array[i,0,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)
                            t = secs[peaks[i]:troughs[i+1]]
                            power = P_in(range_array[i,1,j,k],t,gamma[j,k])
                            for p in power:
                                power_density.append(p)

                        t = secs[troughs[-2]:peaks[-1]]
                        power = P_in(range_array[i+1,0,j,k],t,gamma[j,k])
                        for p in power:
                            power_density.append(p)
                        t = secs[peaks[-1]:troughs[-1]+1]
                        power = P_in(range_array[i+1,1,j,k],t,gamma[j,k])
                        for p in power:
                            power_density.append(p)
                            
                        t_dt = dt[troughs[0]:troughs[-1]+1]
                            
                            


                power_density = np.array(power_density)            
                pd[:len(power_density),j,k] = power_density
                
                t_pd[:len(t_dt),j,k] = t_dt



    return pd,t_pd









def write_pd(pd,pd_time,latitude,longitude,filename):
    '''
    '''
    power_xr = xr.Dataset({
        'pd': xr.DataArray(
            data=pd,
            coords={'series': np.arange(1,pd.shape[0]+1),
                    'lat': latitude,
                    'lon': longitude},
            dims=['series','lat','lon']),
        
        'time': xr.DataArray(
            data=pd_time,
            coords={'series': np.arange(1,pd.shape[0]+1),
                    'lat': latitude,
                    'lon': longitude},
            dims=['series','lat','lon'])})
                                   

    power_xr.to_netcdf(path=f'data/processed/{filename}')
    
    return

    

