%% a0_parameters_setup.m
%  Master parameter file for SESAME HBI posterior surface inversion
%  Adapted from Brunsvik et al. ENAM setup

clear; close all;

%% Paths
results_dir = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/Results_and_Figures/real_data_correlated_results';
datapack_path = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/DATA/sta_datapack.mat';
summary_file = fullfile(results_dir, 'median_vs_summary.mat');
out_dir = fullfile(results_dir, '3d_volume');
func_dir = '/Volumes/AS-Filer/EES/jbrussel/SharedData/SESAME/calc_HBI_posterior_jbr/functions';
fig_dir = fullfile(results_dir, '3d_volume', 'figures');

if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
addpath(func_dir);
addpath(genpath('~/Documents/MATLAB/toolboxes/m_map'));

%% Specify details of this run
ifsave = true;
STAMP = 'correlated_30k';

%% Quality thresholds
desired_chains = 10;
desired_iter = 24000;

%% Load data
load(summary_file);
load(datapack_path);
nsta = length(summary);
sta_lat = [summary.lat]';
sta_lon = [summary.lon]';
sta_names = {summary.name}';

%% Depth grid for surface inversion
zcoarse = summary(1).zcoarse;
vs_vec_pdf = summary(1).vs_vec_pdf;
nz = length(zcoarse);
nvs = length(vs_vec_pdf);
z_vs = zcoarse;

%% Collect station median profiles and PDFs
med_vs_all = zeros(nsta, nz);
pdf_all = zeros(nsta, nz, nvs);
for is = 1:nsta
    med_interp = interp1(summary(is).zinterp, summary(is).med_vs, zcoarse, 'linear', NaN);
    med_vs_all(is,:) = med_interp;
    pdf_all(is,:,:) = summary(is).marginal_pdf;
end

%% Vertical smoothing of station profiles
smooth_nz = 11;
for is = 1:nsta
    med_vs_all(is,:) = movmean(med_vs_all(is,:), smooth_nz, 'omitnan');
end

%% Interface medians and PDFs
med_zmoh_all = [summary.med_zmoh]';
med_zsed_all = [summary.med_zsed]';

%% Regular lat/lon grid
dlat = 0.05;
dlon = 0.05;
lat_pad = 0.1;
lon_pad = 0.1;
lonmin = min(sta_lon) - lon_pad;
lonmax = max(sta_lon) + lon_pad;
latmin = min(sta_lat) - lat_pad;
latmax = max(sta_lat) + lat_pad;
llminmax = [lonmin, lonmax, latmin, latmax];

lat_vec = (latmin : dlat : latmax)';
lon_vec = (lonmin : dlon : lonmax)';
[LON, LAT] = meshgrid(lon_vec, lat_vec);
nlat = length(lat_vec);
nlon = length(lon_vec);

%% Map projection (Lambert) for grid coordinates
m_proj('lambert', 'long', [lonmin-0.5, lonmax+0.5], 'lat', [latmin-0.5, latmax+0.5]);
[stax, stay] = m_ll2xy(sta_lon, sta_lat);
[xgrid, ygrid] = m_ll2xy(LON, LAT);

%% Gaussian interpolation length scale
pair_dist = zeros(nsta);
for i = 1:nsta
    for j = i+1:nsta
        pair_dist(i,j) = deg2km(distance(sta_lat(i), sta_lon(i), sta_lat(j), sta_lon(j)));
        pair_dist(j,i) = pair_dist(i,j);
    end
end
pair_dist(pair_dist == 0) = Inf;
nearest_dist = min(pair_dist, [], 2);
median_spacing = median(nearest_dist);
sigma_gauss = median_spacing * 1.5;

%% Roughness penalty parameters (following Brunsvik et al.)
rough_scale_base = 1e-9;
% rough_scale_params = struct(...
%     'zmoh', 0.0015 * rough_scale_base, ...
%     'zsed', 0.1 * rough_scale_base);

rough_scale_params = struct(...
    'zmoh', 0.00001 * rough_scale_base, ...
    'zsed', 0.001 * rough_scale_base);
%% Inversion parameters
max_inv_iterations = 30;
re_run = false;
vs_clamp = [0.1, 5.5];
zmoh_clamp = [5, 80];
zsed_clamp = [0.001, 5.0];

%% Map plotting
try
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    have_states = true;
catch
    have_states = false;
end

map_lat = [min(sta_lat)-1.0, max(sta_lat)+1.0];
map_lon = [min(sta_lon)-1.0, max(sta_lon)+1.0];
marker_sz = 100;

state_labels = {
    'GA', -83.5, 32.6;
    'SC', -80.9, 33.9;
    'NC', -80.5, 35.5;
    'TN', -85.5, 35.8;
};

%% Output file paths
fresults = fullfile(out_dir, sprintf('compiled_results_%s.mat', STAMP));
fpdfs = fullfile(out_dir, sprintf('compiled_pdfs_%s.mat', STAMP));

fprintf('Setup complete: %d stations, %d depths, grid %dx%d\n', nsta, nz, nlat, nlon);
fprintf('Median spacing: %.1f km, Gaussian sigma: %.1f km\n', median_spacing, sigma_gauss);
fprintf('Results dir: %s\n', results_dir);
fprintf('Output dir: %s\n', out_dir);