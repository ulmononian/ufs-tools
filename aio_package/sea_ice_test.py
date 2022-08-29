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
import cmocean.cm as cmo


land_50m = cft.NaturalEarthFeature('physical', 'land', '50m',
                                   edgecolor='black', facecolor='papayawhip', linewidth=0.5)

def get_circle():
    """
    Compute a circle in axes coordinates, which we can use as a boundary
    for the map. We can pan/zoom as much as we like - the boundary will be
    permanently circular.
    """
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T

    return mpath.Path(verts * radius + center)

def draw_map_subplots(ax):
    projection = "PolarNorth"
    ax.set_global()
    ax.add_feature(cartopy.feature.LAND,zorder=1)
    ax.add_feature(cartopy.feature.COASTLINE,zorder=1)
    ax.set_aspect('auto',adjustable=None)
    trans = ccrs.PlateCarree()
    if projection == 'PolarNorth':
        ax.set_extent([-180,180,60,90],crs=ccrs.PlateCarree())
        kw=dict(central_latitude=90,central_longitude=-45,true_scale_latitude=70)
    if projection == 'PolarSouth':
        ax.set_extent([-180,180,-90,-60],crs=ccrs.PlateCarree())
        kw=dict(central_latitude=-90,central_longitude=0,true_scale_latitude=-70)
    if projection == 'PolarSouth' or projection == 'PolarNorth':
        ax.gridlines()
        circle_value = get_circle()
        ax.set_boundary(circle_value, transform=ax.transAxes)
        trans = ccrs.Stereographic(**kw)
    return ax, trans

ocn_file = "ocn_2021_03_22_09.nc"
print(ocn_file)
#ocn = xr.open_dataset(ocn_file,decode_times=False)
#yh = ocn[yh][:]; xh = ocn[xh][:]
ocn = nc.Dataset(ocn_file)
yh = ocn.variables['yh'][:]
xh = ocn.variables['xh'][:]
geolon = ocn.variables['geolon'][:]
geolat = ocn.variables['geolat'][:]


ice_file = nc.Dataset('iced.2021-03-23-43200.nc','r')
field_ =  np.ma.masked_invalid(ice_file['aicen'][:])
#field_ = ice_file.variables['aicen'][:]
#field1 = np.ma.masked_invalid(field_['aicen'])
field1 = (field_[0,:,:]+field_[1,:,:]+field_[2,:,:]+field_[3,:,:]+field_[4,:,:])/5.

print(np.amax(field1))
print(np.amin(field1))
with np.printoptions(threshold=np.inf):
    print(field1)

min = np.amin(field1)
max = np.amax(field1)
levels= np.arange(np.amin(field1),np.amax(field1), 0.01)
# min_  = levels[0]
# max_  = levels[-1] + field_range[2]
#
#fig = plt.figure()
#ax = plt.axes(projection=ccrs.SouthPolarStereo())
#ax.add_feature(cartopy.feature.LAND)
#ax.add_feature(cartopy.feature.COASTLINE)
#ax.set_aspect('auto',adjustable=None)
#ax.set_extent([-180,180,60,90],crs=ccrs.PlateCarree())
#kw=dict(central_latitude=90,central_longitude=-45,true_scale_latitude=70)
#ax.set_extent([-180,180,-90,-60],crs=ccrs.PlateCarree())
#kw=dict(central_latitude=-90,central_longitude=0,true_scale_latitude=-70)
#theta = np.linspace(0, 2*np.pi, 100)
#center, radius = [0.5, 0.5], 0.5
#verts = np.vstack([np.sin(theta), np.cos(theta)]).T
#circle = mpath.Path(verts * radius + center)
#ax.set_boundary(circle, transform=ax.transAxes)
#trans = ccrs.Stereographic(**kw)
#cs = plt.contourf(geolon,geolat, field1, cmap=cmo.ice, transform=ccrs.PlateCarree())
#cb = fig.colorbar(cs, orientation='vertical')
projection=ccrs.SouthPolarStereo()

fig = plt.figure(figsize=(10, 9))
ax = plt.subplot(1, 1, 1, projection=projection)

ax.set_extent([-180,180,-90,-60], crs=ccrs.PlateCarree())
ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])
ax.coastlines(resolution='50m')

# Compute a circle in axes coordinates, which we can use as a boundary
# for the map. We can pan/zoom as much as we like - the boundary will be
# permanently circular.
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)

p1 = ax.pcolormesh(geolon, geolat, field1, vmin=min, vmax=max, cmap=cmo.ice, transform=ccrs.PlateCarree())

ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
cb = plt.colorbar(p1, cax=ax_cb, orientation='vertical')
#cb.ax.set_ylabel('Sea Ice Concentration');
plt.savefig('test2.png')
