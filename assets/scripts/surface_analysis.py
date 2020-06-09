import titanlib
import netCDF4
import numpy as np

# Load observations
with netCDF4.Dataset('obs.nc', 'r') as file:
    obs_lats = file.variables['latitude'][:]
    obs_lons = file.variables['longitude'][:]
    obs_elevs = file.variables['altitude'][:]
    obs = file.variables['air_temperature_2m'][:, 0]

# SCT settings
inner_radius = 50000
outer_radius = 100000
num_min = 5
num_max = 100
num_iterations = 1
num_min_prof = 20
dzmin = 100
dhmin = 10000
dz = 200
t2pos = np.full([len(obs)], 4)
t2neg = np.full([len(obs)], 4)
eps2 = np.full([len(obs)], 0.5)
print(t2pos.shape)

# Run the SCT
flags, sct, rep = titanlib.sct(obs_lats, obs_lons, obs_elevs, obs, num_min, num_max, inner_radius,
        outer_radius, num_iterations, num_min_prof, dzmin, dhmin , dz, t2pos, t2neg, eps2)

index_valid_obs = np.where(flags == 0)[0]
index_invalid_obs = np.where(flags != 0)[0]
print(index_invalid_obs)

# Plot the quality control results
import matplotlib.pylab as mpl
mpl.scatter(obs_lons[index_valid_obs], obs_lats[index_valid_obs],
        c=obs[index_valid_obs] - 273.15, s=50, linewidths=0, cmap="RdBu_r", vmin=5, vmax=20)
mpl.scatter(obs_lons[index_invalid_obs], obs_lats[index_invalid_obs],
        c=obs[index_invalid_obs] - 273.15, s=50, edgecolors="r", linewidths=1.5, cmap="RdBu_r",
        vmin=5, vmax=20)
mpl.xlim(0, 35)
mpl.ylim(55, 75)
mpl.gca().set_aspect(2)
cb = mpl.colorbar()
cb.set_label(r"Temperature ($\degree C$)")
mpl.savefig("titan_sct.png", bbox_inches='tight', dpi=125)
mpl.clf()

# Read NWP background from file
import gridpp
import numpy as np

with netCDF4.Dataset('analysis.nc', 'r') as file:
    index_control = 0
    blats_2500m = file.variables['latitude'][:]
    blons_2500m = file.variables['longitude'][:]
    belevs_2500m = file.variables['surface_geopotential'][0, 0, index_control, :, :] / 9.81
    bgrid_2500m = gridpp.Grid(blats_2500m, blons_2500m, belevs_2500m)
    background_2500m = file.variables['air_temperature_2m'][0, 0, index_control, :, :]

with netCDF4.Dataset('template.nc', 'r') as file:
    blats = file.variables['latitude'][:]
    blons = file.variables['longitude'][:]
    belevs = file.variables['altitude'][:]
    bgrid = gridpp.Grid(blats, blons, belevs)

# Downscaling
gradient = -0.0065
background = gridpp.simple_gradient(bgrid_2500m, bgrid, background_2500m, gradient)

points = gridpp.Points(obs_lats[index_valid_obs], obs_lons[index_valid_obs], obs_elevs[index_valid_obs])
pbackground = gridpp.bilinear(bgrid, points, background)

variance_ratios = np.full(points.size(), 0.1)

h = 100000
v = 200
structure = gridpp.BarnesStructure(h, v)

max_points = 50

analysis = gridpp.optimal_interpolation(bgrid, background, points,
        obs[index_valid_obs], variance_ratios, pbackground, structure, max_points)


# Plotting the increments
diff = analysis - background
mpl.pcolormesh(blons, blats, diff, cmap="RdBu_r", vmin=-2, vmax=2)
mpl.xlim(0, 35)
mpl.ylim(55, 75)
mpl.gca().set_aspect(2)
cb = mpl.colorbar()
cb.set_label(r"Increment ($\degree C$)")
mpl.savefig("analysis_increment.png", bbox_inches='tight', dpi=125)
mpl.clf()

# Ensemble mode
with netCDF4.Dataset('analysis.nc', 'r') as file:
    background_ens_2500m = file.variables['air_temperature_2m'][0, 0, :, :, :]
    background_ens_2500m = np.moveaxis(background_ens_2500m, 0, 2)

num_members = background_ens_2500m.shape[2]
pbackground_ens = np.zeros([points.size(), num_members])
background_ens = np.zeros([bgrid.size()[0], bgrid.size()[1], num_members])
for e in range(num_members):
    background_ens[:, :, e] = gridpp.simple_gradient(bgrid_2500m, bgrid, background_ens_2500m[:, :, e], gradient)
    pbackground_ens[:, e] = gridpp.bilinear(bgrid, points, background_ens[:, :, e])
print("Done downscaling")

psigmas = np.full(points.size(), 0.5)

analysis_ens = gridpp.optimal_interpolation_ensi(bgrid, background_ens, points,
        obs[index_valid_obs], psigmas, pbackground_ens, structure, max_points)
diff = (analysis_ens - background_ens)[:, :, index_control]

mpl.pcolormesh(blons, blats, diff, cmap="RdBu_r", vmin=-2, vmax=2)
mpl.xlim(0, 35)
mpl.ylim(55, 75)
mpl.gca().set_aspect(2)
cb = mpl.colorbar()
cb.set_label(r"Increment ($\degree C$)")
mpl.savefig("analysis_increment_ens.png", bbox_inches='tight', dpi=125)
