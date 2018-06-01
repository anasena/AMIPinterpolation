#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:55:16 2018

@author: athomese

Calculating the mid-monthly SST from the monthly mean SST.
Based on https://pcmdi.llnl.gov/mips/amip/details/#Appendix_1
"""
import numpy as np
from netCDF4 import Dataset

def getting_dimensions(filein):
    lats=dataset.variables["lat"][:]
    lons=dataset.variables["lon"][:]
    time=dataset.variables["time"][:]
    date=dataset.variables["date"][:]
    datesec=dataset.variables["datesec"][:]
    icecov=dataset.variables["ice_cov"][:]
    sst=dataset.variables["SST_cpl"][:]
    # In the case climatology case, I only have one year. So, I can replicate to minimize the boundary condition problems:
    if time.size==12:
        time=np.arange(0,24,1)
        date=np.concatenate((date,date))
        datesec=np.concatenate((datesec,datesec))
        icecov=np.concatenate((icecov,icecov),axis=0)
        sst=np.concatenate((sst,sst),axis=0)
    return(lats,lons,time,date, datesec, icecov, sst)


def print_netcdf(variables, var_names, filename, filein):
    print('\n>>> Using netcdf function!')
    lats,lons,time,date, datesec, icecov, sst_original=getting_dimensions(filein)
    
    for i in range(0, len(var_names)):
        var_name=var_names[i]
        var=variables[i]
        if i==0:
            ncfile = Dataset(filename,'w',format='NETCDF4')
        #    ncfile.description = ' ' \
        #                        ' ' 
            
            ### Dimensions
            ncfile.createDimension('time',var.shape[0])    
            ncfile.createDimension('lat',var.shape[1])
            ncfile.createDimension('lon',var.shape[2])
            
            ### Variables
            latitude = ncfile.createVariable('lat','f4',('lat',))
            longitude = ncfile.createVariable('lon','f4',('lon',))
            times = ncfile.createVariable('time','f4',('time',))
            dates = ncfile.createVariable('date','f4',('time',))
            datesecs = ncfile.createVariable('datesec','f4',('time',))
            icecovs=ncfile.createVariable('ice_cov_prediddle','f4',('time','lat','lon'))
            sst=ncfile.createVariable('SST_cpl_prediddle','f4',('time','lat','lon'))
            latitude[:] = lats
            longitude[:] = lons
            times[:]=time
            dates[:]=date
            datesecs[:]=datesec
            icecovs[:]=icecov
            sst[:]=sst_original          
    # Here goes your dataset from the correction:
        varns = ncfile.createVariable(var_name,'f4',('time','lat','lon'))
    ### Data
        varns[:] = var
    ### Units
#    varns.units = '%'
#    ncfile.title = ''
#    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
#    ncfile.source = ''
#    ncfile.references = 'SSMIS (DMSP F18)]'
    
    ncfile.close()
    print('*Completed: Created netCDF4 File!')


def interpolating_matrix(S):
    # input: S = observed mean monthly SST in one grid point
    # output: A=The matrix to defining the equation system
    S=S.squeeze()
    A=np.zeros((S.size,S.size))
    # boundary conditions:
    # in the cases you can ignore the beginning of your timeseries:
    # The first temperature (T[0]) will be equal to the monthly mean, as well as the final condition:
    A[0,0]=1
    A[-1,-1]=1

    for i in range (1,S.size-1):
        A[i,i-1:i+1+1]=[1/8,3/4,1/8]
    return A



# insert the name of the files here:
directory_in='/surtsey/ypeings/'
directory_out='/surtsey/athomese/'
var_names=["SST_cpl", "ice_cov"]
filenames=['SST_hadisst_clim_1979-2008.nc', 'SST_hadisst_IODp_cesm1.2.nc', 'SST_hadisst_IODn_cesm1.2.nc']

for k in range (0, len(filenames)):
    filename=filenames[k]
    dataset=Dataset(directory_in + filename)
    fileout=directory_out+filename[0:-3]+'_AMIP_correction.nc'

    T_vars=[]
    for z in range (0, len(var_names)):
        var_name=var_names[z];        
    
        S=dataset.variables[var_name][:]
    
        if filename=='SST_hadisst_clim_1979-2008.nc':
            S=np.concatenate((S,S), axis=0)
        # for each gridpoint:
        # input: S = observed mean monthly SST in one grid point
        # output: T=mid-month SST in that grid point
    
        T=np.empty(S.shape)
        for i in range(0,S.shape[1]):
            for j in range(0,S.shape[2]):
                A=interpolating_matrix(S[:,i,j])
                # I'll determine the new midpoint temperature
                T[:,i,j]=np.linalg.solve(A,S[:,i,j])            
        T_vars.append(T)
    print_netcdf(T_vars, var_names, fileout, dataset)
    dataset.close()
