%% Gaussian-weighted interface surfaces (Moho and Sediment)
%  Use instead of fminunc for interfaces since derived PDFs are too broad

run('a0_parameters_setup.m');

sfs = load(fullfile(out_dir, 'surface_out_example.mat'));
xgrid_inv = sfs.xgrid;
ygrid_inv = sfs.ygrid;
longrid = sfs.longrid;
latgrid = sfs.latgrid;

%% Compute distance from each grid node to each station
dist_mat = zeros(size(longrid,1), size(longrid,2), nsta);
for is = 1:nsta
    dist_mat(:,:,is) = deg2km(distance(latgrid, longrid, sta_lat(is), sta_lon(is)));
end
max_dist = 3 * sigma_gauss;
weights = exp(-dist_mat.^2 / (2*sigma_gauss^2));
weights(dist_mat > max_dist) = 0;

%% Gaussian-weighted Moho surface
valid_moh = ~isnan(med_zmoh_all);
w_moh = weights(:,:,valid_moh);
wsum_moh = sum(w_moh, 3);
zmoh_weighted = zeros(size(longrid));
idx_moh = find(valid_moh);
for k = 1:length(idx_moh)
    zmoh_weighted = zmoh_weighted + w_moh(:,:,k) * med_zmoh_all(idx_moh(k));
end
good = wsum_moh > 1e-10;
zmoh_gauss = NaN(size(longrid));
zmoh_gauss(good) = zmoh_weighted(good) ./ wsum_moh(good);

%% Gaussian-weighted Sediment surface
valid_sed = ~isnan(med_zsed_all);
w_sed = weights(:,:,valid_sed);
wsum_sed = sum(w_sed, 3);
zsed_weighted = zeros(size(longrid));
idx_sed = find(valid_sed);
for k = 1:length(idx_sed)
    zsed_weighted = zsed_weighted + w_sed(:,:,k) * med_zsed_all(idx_sed(k));
end
good_sed = wsum_sed > 1e-10;
zsed_gauss = NaN(size(longrid));
zsed_gauss(good_sed) = zsed_weighted(good_sed) ./ wsum_sed(good_sed);

%% Save to interface files
mgrid_out = zmoh_gauss;
save(fullfile(out_dir, 'zmoh', 'surface_values_V1.mat'), 'mgrid_out');
mgrid_out = zsed_gauss;
save(fullfile(out_dir, 'zsed', 'surface_values_V1.mat'), 'mgrid_out');

%% Update 3D volume
vol_data = load(fullfile(out_dir, 'vs_3d_volume.mat'));
vol_data.zmoh_inv = zmoh_gauss;
vol_data.zsed_inv = zsed_gauss;
save(fullfile(out_dir, 'vs_3d_volume.mat'), '-struct', 'vol_data', '-v7.3');

fprintf('zmoh Gaussian: min=%.1f max=%.1f range=%.1f\n', ...
    nanmin(zmoh_gauss(:)), nanmax(zmoh_gauss(:)), nanmax(zmoh_gauss(:))-nanmin(zmoh_gauss(:)));
fprintf('zsed Gaussian: min=%.3f max=%.3f range=%.3f\n', ...
    nanmin(zsed_gauss(:)), nanmax(zsed_gauss(:)), nanmax(zsed_gauss(:))-nanmin(zsed_gauss(:)));
fprintf('Updated interface files and 3D volume\n');