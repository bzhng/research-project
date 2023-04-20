### by observation, the stream power time series is rather erratic, thus at this point in time, the peaks have not been found


import os
from time import time
import multiprocessing
import netCDF4 as nc
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas




epoch = np.datetime64('1970-01-01T00:00:00')
T = ((12*60+25)*60)/2 


def peak_trough_idx(height,t,latitude,longitude):
    '''
    indices of the peaks and troughs
    
    '''
#     count = 0
    lat_error = []
    lon_error = []
    
    
    secs = (t-epoch)/1e9 # dtype is timedelta[s] not timedelta[ns]
    secs = secs.astype('float64')
    
    step = 25
    istep = 20
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