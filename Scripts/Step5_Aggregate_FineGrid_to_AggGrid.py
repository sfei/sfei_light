# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 15:55:32 2019
Take the interpolated full-resolution light field and perform a spatial aggregation with the aggregated grid
the module dpp can be found on hpc:/hpcvol1/zhenlin/Ztoolbox
@author: zhenlinz
"""

#import dwaq.PostProcessing as dpp
import geopandas as gpd
from stompy import utils
from stompy.grid import unstructured_grid
import numpy as np
import logging
import xarray as xr
import netCDF4
import os, glob


dir_grid = "../Grid"
file_finegrid = os.path.join(dir_grid,'wy2013c_waqgeom.nc')
file_agggrid = os.path.join(dir_grid,'flowgeom141.nc')

#gridf = dpp.dwaqGrid(file_finegrid)
gridf = unstructured_grid.UnstructuredGrid.read_dfm(file_finegrid)
grida = unstructured_grid.UnstructuredGrid.read_dfm(file_agggrid)

gpdf = gridf.write_cells_geopandas()
gpda = grida.write_cells_geopandas()

mapping = gpd.sjoin(gpdf, gpda, how="left", predicate='within')

# The following should not happen but just a quick check
#if np.any(np.isnan(mapping.index_right.values)):
#    contiguous=False
#    indnan = np.where(np.isnan(mapping.index_right.values))[0]
#    logging.warning("unclassfied cells identified: correction needed at {}".format(indnan))
    
# %% Some really messed up workspace clearing to handle xarray bugs
# RH: hoping this is no longer necessary   
#for name in dir():
#    if not (name.startswith('file_output')) | (name.startswith('_')) | (name.startswith('gpd')) | (name.startswith('file')) | (name.startswith('mapping')):
#        del globals()[name]
        
import xarray as xr
import numpy as np
import os

dir_input = '../Data_DELWAQ_Inputfiles_RH'
dir_output = '../Data_DELWAQ_Agg_Inputfiles_RH'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)
    
input_names = glob.glob(os.path.join(dir_input,'Kd_PropShift_forDELWAQ_*.nc'))
for file_input in input_names:
    print("----- %s -----"%os.path.basename(file_input))
    
    file_output = os.path.join(dir_output,os.path.basename(file_input)[:-3]+'_AggGrid_141.nc')

    light = xr.open_dataset(file_input)

    # original code loaded the entire Kd dataset at once.
    # that could be 30+GB.
    # Better to iterate over time

    if 0:
        light_full = light['Kd'].values # Hmm - this is potentially very large!
        light_avg = []
        for i in np.arange(len(gpda)): # for each aggregated polygon
            poly_cells = np.nonzero(mapping.index_right==i)[0] # the included high resolution grid cells
            light_i = light_full[:,poly_cells].mean(axis=1)
            light_avg.append(light_i)


        time = light.time.values
        nFlowElem = gpda.index.values
        light_avg = np.asarray(light_avg).T
        light_new = xr.DataArray(light_avg,coords=[time,nFlowElem],dims=['time','nFlowElem'])
        light_new = light_new.to_dataset(name='kd')

        light_new.to_netcdf(file_output)
    else:
        ds_new=xr.Dataset()
        ds_new['time']=light.time
        # generally nice to keep the grid geometry nearby
        grida.write_xarray(ds=ds_new,face_dimension='nFlowElem',edge_dimension='nNetLink',node_dimension='nNetNode')
        ds_new.to_netcdf(file_output)
        ds_new.close()
        del ds_new
        nc=netCDF4.Dataset(file_output,mode='a')
        kd_var=nc.createVariable('Kd',"f8",['time','nFlowElem'])

        #light_full = light['Kd'].values # Hmm - this is potentially very large!
        #light_avg = []

        poly_cells={} # agg polygon => array of fine grid cell indices
        for i in np.arange(len(gpda)): # for each aggregated polygon
            poly_cells[i] = np.nonzero(mapping.index_right.values==i)[0] # the included high resolution grid cells

        # Turns out it's pretty slow to load a single time index at a time. Not much extra cost
        # to load 100 (or more). So go through, reading 100 times, aggregate, write out.
        stride=100
        def slices(N,stride=100):
            i=0
            while i<N:
                i_next=min(i+stride,N)
                yield slice(i,i_next)
                i=i_next
                
        for ti in slices(light.dims['time']):
            print(ti)
            # area-weighted would be better, but it's not going to make a big difference.
            # this could be streamlined, but probably the cost is all in in reading the data off disk.
            #light_snap=light['Kd'].isel(time=ti).values

            light_snap=light['Kd'].isel(time=ti).values # slices, so [times,cells]
            
            light_snap_agg=[light_snap[...,poly_cells[i]].mean(axis=-1)
                            for i in np.arange(len(gpda))]
            # => [polygons, times]
            
            kd_var[ti,:]=np.array(light_snap_agg).T
        nc.close()


