run('a0_parameters_setup.m');

% List output directories
d = dir(out_dir);
d = d([d.isdir] & ~ismember({d.name},{'.','..'}));
fprintf('Output directories:\n');
for i = 1:length(d)
    has_result = exist(fullfile(out_dir, d(i).name, 'surface_values_V1.mat'), 'file');
    if has_result
        fprintf('  %s: done\n', d(i).name);
    else
        fprintf('  %s: missing\n', d(i).name);
    end
end

% Quick check: load and plot a few depth slices
check_depths = [0.5, 5, 30, 50];
figure('Position',[50 50 1200 800],'Color','w');
for ip = 1:length(check_depths)
    this_inv = sprintf('vs%.1f', check_depths(ip));
    f = fullfile(out_dir, this_inv, 'surface_values_V1.mat');
    if ~exist(f,'file'), fprintf('Missing %s\n', this_inv); continue; end
    mgrid_out = load(f).mgrid_out;
    subplot(2,2,ip); hold on;
    contourf(longrid, latgrid, mgrid_out, 15, 'LineStyle','none');
    plot(sta_lon, sta_lat, 'k^', 'MarkerFaceColor','w', 'MarkerSize',5);
    colorbar; colormap(turbo);
    title(sprintf('Vs at %.1f km', check_depths(ip)));
    axis equal tight;
end
sgtitle('Surface Inversion Output Check');
exportgraphics(gcf, fullfile(out_dir, 'diagnostic_check.png'), 'Resolution', 300);
fprintf('Saved diagnostic\n');



run('a0_parameters_setup.m');

sfs = load(fullfile(out_dir, 'surface_out_example.mat'));
longrid = sfs.longrid;
latgrid = sfs.latgrid;

check_depths = [0.5, 5, 30, 50];
figure('Position',[50 50 1200 800],'Color','w');
for ip = 1:length(check_depths)
    this_inv = sprintf('vs%.1f', check_depths(ip));
    f = fullfile(out_dir, this_inv, 'surface_values_V1.mat');
    if ~exist(f,'file'), fprintf('Missing %s\n', this_inv); continue; end
    mgrid_out = load(f).mgrid_out;
    subplot(2,2,ip); hold on;
    contourf(longrid, latgrid, mgrid_out, 15, 'LineStyle','none');
    plot(sta_lon, sta_lat, 'k^', 'MarkerFaceColor','w', 'MarkerSize',5);
    colorbar; colormap(turbo);
    title(sprintf('Vs at %.1f km', check_depths(ip)));
    axis equal tight;
end
sgtitle('Surface Inversion Output Check');
exportgraphics(gcf, fullfile(out_dir, 'diagnostic_check.png'), 'Resolution', 300);
fprintf('Saved diagnostic\n');