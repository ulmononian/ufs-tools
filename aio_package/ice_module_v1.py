import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cft
import cmocean
import cmocean.cm as cm
import json


def get_circle():
    """
    Compute a circle in axes coordinates, which we can use as a boundary
    for the map. We can pan/zoom as much as we like - the boundary will be
    permanently circular.
    """
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle_value = mpath.Path(verts * radius + center)
    return circle_value

def draw_map_subplots(ax,d):
    projection = d["projection"]
    trans = ccrs.PlateCarree()
    if projection == 'PolarNorth':
        ax.set_extent([-180,180,60,90],crs=ccrs.PlateCarree())
        kw=dict(central_latitude=90,central_longitude=-45,true_scale_latitude=70)
    if projection == 'PolarSouth':
        ax.set_extent([-180,180,-90,-60],crs=ccrs.PlateCarree())
        kw=dict(central_latitude=-90,central_longitude=0,true_scale_latitude=-70)
    if projection == 'PolarSouth' or projection == 'PolarNorth':
        #ax.gridlines()
        circle_value = get_circle()
        ax.set_boundary(circle_value, transform=ax.transAxes)
        trans = ccrs.Stereographic(**kw)
    return ax, trans

def set_projection_cmap(projection, cmap_in):
    if projection == 'PolarNorth': crs=ccrs.NorthPolarStereo()
    if projection == 'PolarSouth': crs=ccrs.SouthPolarStereo()
    if projection == 'Robinson'  : crs=ccrs.Robinson(central_longitude=-120)

    if cmap_in == "plt.cm.jet": cmap_=cm.jet
    #if cmap_in == "plt.cm.jet": cmap_=cmap_map(lambda x:x*.85,matplotlib.cm.jet)
    if cmap_in == "plt.cm.coolwarm": cmap_=cm.coolwarm
    if cmap_in == "plt.cm.seismic": cmap_=cm.seismic
    if cmap_in == "cmo.RdBu_r": cmap_=cm.RdBu_r
    if cmap_in == "cmo.thermal": cmap_=cm.thermal
    if cmap_in == "cmo.ice": cmap_=cm.ice
    return crs, cmap_

def draw_colorbar(ax, cs, d, fig):
    cbar_label = d["cbar_label"]
    #cbar_label_font = d["cbar_label_font"]
    cax,kw1=matplotlib.colorbar.make_axes([ax],location='bottom',shrink=0.6, pad=0.05, fraction=0.1)
    cbar=fig.colorbar(cs,cax=cax,**kw1)
    #cbar.ax.tick_params(labelsize=cbar_label_font)
    cbar.set_label(cbar_label,size=6)

def make_ice_map(json_filename):
    with open(json_filename) as files:
      data = json.load(files)
    d = data[0]
    ice_file  = d["ice_file"]
    grid   = d["grid"]
    figsize= d["figsize"]
    fields = d["fields"]
    projection = d["projection"]
    cmap_in    = d["cmap"]
    title_ = d["title"]
    title_font_ = d["title_font"]
    filename=d["output_fig"]

    # get grid info from MOM6 geolon,geolat
    ocn = nc.Dataset(grid)
    geolon = ocn.variables['geolon'][:]
    geolat = ocn.variables['geolat'][:]

    # load CICE file and get field information
    ice_in = nc.Dataset(ice_file[0],'r')
    field_ = ice_in.variables[fields[0]]
    field1 = (field_[0,:,:]+field_[1,:,:]+field_[2,:,:]+field_[3,:,:]+field_[4,:,:])/5.

    field1_min = np.amin(field1)
    field1_max = np.amax(field1)

    fig = plt.figure(figsize=figsize)

    crs, cmap_ = set_projection_cmap(projection, cmap_in)
    ax = plt.subplot(1,1,1, projection=crs)
    ax,trans = draw_map_subplots(ax,d)
    land_50m = cft.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='papayawhip', linewidth=0.5)
    ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
    ax.coastlines(resolution='50m')

    cs = ax.pcolormesh(geolon, geolat, field1, vmin=field1_min, vmax=field1_max, cmap=cmap_, transform=ccrs.PlateCarree())
    cb = draw_colorbar(ax,cs,d,fig)
    ax.set_title(title_,fontsize=title_font_)
    plt.savefig('test5.png')

