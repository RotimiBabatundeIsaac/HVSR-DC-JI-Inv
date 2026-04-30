%% a4_collate_3d_volume.m
%  Assemble individual surface inversion outputs into a 3D volume
%  Then slice cross-sections through it

run('a0_parameters_setup.m');

sfs = load(fullfile(out_dir, 'surface_out_example.mat'));
longrid = sfs.longrid;
latgrid = sfs.latgrid;
xgrid = sfs.xgrid;
ygrid = sfs.ygrid;

version_surf = 1;

%% Load all depth slices into 3D volume
vol_inv = NaN(size(longrid,1), size(longrid,2), nz);
for iz = 1:nz
    this_inv = sprintf('vs%.1f', zcoarse(iz));
    f = fullfile(out_dir, this_inv, sprintf('surface_values_V%d.mat', version_surf));
    if exist(f, 'file')
        vol_inv(:,:,iz) = load(f).mgrid_out;
    else
        fprintf('Missing %s\n', this_inv);
    end
end
fprintf('Loaded %d depth slices\n', sum(~isnan(vol_inv(1,1,:))));

%% Load interface surfaces
f_zmoh = fullfile(out_dir, 'zmoh', sprintf('surface_values_V%d.mat', version_surf));
f_zsed = fullfile(out_dir, 'zsed', sprintf('surface_values_V%d.mat', version_surf));
zmoh_inv = load(f_zmoh).mgrid_out;
zsed_inv = load(f_zsed).mgrid_out;

%% Save assembled volume
save(fullfile(out_dir, 'vs_3d_volume.mat'), ...
    'vol_inv', 'zmoh_inv', 'zsed_inv', ...
    'longrid', 'latgrid', 'xgrid', 'ygrid', ...
    'zcoarse', 'sta_lat', 'sta_lon', ...
    '-v7.3');
fprintf('Saved 3D volume\n');

%% Diagnostic: check value ranges
for iz_check = [1, 10, 60, 100]
    slab = vol_inv(:,:,iz_check);
    vals = slab(~isnan(slab));
    fprintf('z=%.1f km: min=%.2f max=%.2f\n', zcoarse(iz_check), min(vals), max(vals));
end