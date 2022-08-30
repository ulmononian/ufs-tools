import argparse
import numpy as np
import cmocean
import cmocean.cm as cmo
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.path as mpath
import cartopy
import cartopy.crs as ccrs
import xarray as xr
import datetime
import json
#from matplotlib.colors import TwoSlopeNorm
import netCDF4


def make_ocn_map(json_filename):
  #parser = argparse.ArgumentParser()
  #parser.add_argument('-json', nargs='*', required=True, help='tmp.json')
  #args = parser.parse_args()
  #filename= str(args.json[0])
  #print(filename)

  # read in json file to plot ------------------
  with open(json_filename) as files:
      data = json.load(files)
  d = data[0]
  files  = d["files"]
  grid   = d["grid"]
  level  = d["field_range"]
  figsize= d["figsize"]
  fields = d["fields"]
  projection = d["projection"]
  title_ = d["title"]
  filename=d["output_fig"]

  nc= xr.open_dataset(files[0],decode_times=False); lat= nc['yh']; lon= nc['xh']
  [lon,lat]= np.meshgrid(lon,lat)

  file2read = netCDF4.Dataset(files[0],'r')
  field0_ = file2read.variables[fields[0]]
  field0 = field0_[0,:,:]

  if projection == 'Robinson':
    crs=ccrs.Robinson(central_longitude=-120)
  if projection == 'latlon':
    crs=ccrs.PlateCarree(central_longitude=-120)
  fig = plt.figure()
  ax  = plt.subplot(111, projection=crs)
  ax.coastlines(resolution='110m')
  #ax.gridlines()

  print(np.mean(field0),np.std(field0))
  print(np.max(field0),np.min(field0))
  title=title_
  stats='mean:'+str(np.mean(field0))+'\n std:'+str(np.std(field0))

  csst = plt.contourf(lon, lat, field0, levels=np.arange(-4,34,2), cmap=cmo.balance, transform=ccrs.PlateCarree())
  cd0 = fig.colorbar(csst, ax=ax, orientation='horizontal', shrink=0.8, pad=0.03, aspect=50)#, extend="both")
  cd0.ax.locator_params(nbins=10)

  ax.set_title(title+'\n'+stats, fontsize=10)

  #ax.text(-180, -90, stats, ha='left', va='bottom', fontsize=9)
  fig.savefig(filename,dpi=150,facecolor='w',edgecolor='w',transparent=False)

  plt.show()

