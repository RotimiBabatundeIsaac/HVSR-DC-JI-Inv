%% a4_plot_cross_sections.m
%  Slice cross-sections through the 3D Vs volume from surface inversion

run('a0_parameters_setup.m');

vol_data = load(fullfile(out_dir, 'vs_3d_volume.mat'));
vol_inv = vol_data.vol_inv;
zmoh_inv = vol_data.zmoh_inv;
zsed_inv = vol_data.zsed_inv;
longrid = vol_data.longrid;
latgrid = vol_data.latgrid;

basedir = fullfile(out_dir, 'cross_sections');
if ~exist(basedir,'dir'), mkdir(basedir); end

lines = {'D','E','W'};
zoom_depths = [0.5, 5, 70];
zoom_tags = {'500m','5km','70km'};

%% Build 3D coordinate arrays for griddata
nz_vol = length(zcoarse);
lon3d = zeros(size(longrid,1), size(longrid,2), nz_vol);
lat3d = zeros(size(longrid,1), size(longrid,2), nz_vol);
z3d   = zeros(size(longrid,1), size(longrid,2), nz_vol);
for iz = 1:nz_vol
    lon3d(:,:,iz) = longrid;
    lat3d(:,:,iz) = latgrid;
    z3d(:,:,iz) = zcoarse(iz);
end

%% Contour settings (Brennan style)
step_contour = 0.025;
v_contours = (0:step_contour:10) + 0.5*step_contour;
n_colors = 25;

for il = 1:length(lines)
    line_id = lines{il};
    linedir = fullfile(basedir, sprintf('line_%s', line_id));
    if ~exist(linedir,'dir'), mkdir(linedir); end

    %% Identify and sort stations (descending = north to south)
    idx = []; nums = [];
    for is = 1:length(summary)
        nm = summary(is).name;
        if ~isempty(nm) && nm(1) == line_id
            num_val = str2double(nm(2:end));
            if ~isnan(num_val)
                idx(end+1) = is;
                nums(end+1) = num_val;
            end
        end
    end
    [nums, sort_order] = sort(nums, 'descend');
    idx = idx(sort_order);
    nsta_line = length(idx);
    if nsta_line < 2, continue; end
    fprintf('Line %s: %d stations\n', line_id, nsta_line);

    slat = arrayfun(@(i) summary(i).lat, idx)';
    slon = arrayfun(@(i) summary(i).lon, idx)';
    names = arrayfun(@(i) summary(i).name, idx, 'UniformOutput', false)';

    %% Great-circle arc from first to last station
    Q1 = [slat(1), slon(1)];
    Q2 = [slat(end), slon(end)];
    [profd, profaz] = distance(Q1(1), Q1(2), Q2(1), Q2(2));
    nxy = 200;
    gcarc = linspace(0, profd, nxy)';
    d2km = 2*pi*6371/360;
    dist_arc = d2km * gcarc;
    [lat_line, lon_line] = reckon(Q1(1), Q1(2), gcarc, profaz);

    %% Project stations onto section
    sta_dist = zeros(nsta_line, 1);
    for is = 1:nsta_line
        dists = distance(slat(is), slon(is), lat_line, lon_line);
        [~, closest] = min(dists);
        sta_dist(is) = dist_arc(closest);
    end

    %% Slice interfaces
    zmoh_sect = griddata(longrid, latgrid, zmoh_inv, lon_line, lat_line);
    zsed_sect = griddata(longrid, latgrid, zsed_inv, lon_line, lat_line);

    for iz_plot = 1:3
        zmax_plot = zoom_depths(iz_plot);
        z_sub = zcoarse(zcoarse <= zmax_plot);
        nz_sub = length(z_sub);
        if nz_sub < 2, continue; end

        %% Slice volume along section
        [lonmesh, zmesh] = ndgrid(lon_line, z_sub);
        [latmesh, ~]     = ndgrid(lat_line, z_sub);
        [gcmesh, ~]      = ndgrid(gcarc, z_sub);

        Vg = griddata(lon3d, lat3d, z3d, vol_inv, lonmesh, latmesh, zmesh);
        dist_sect = gcmesh * d2km;

        vs_vals = Vg(~isnan(Vg));
        if isempty(vs_vals), continue; end

        %% Mask crust vs mantle
        zmoh_mat = repmat(zmoh_sect(:), 1, nz_sub);
        is_crust  = zmesh <= zmoh_mat;
        is_mantle = zmesh >  zmoh_mat;

        Vg_crust = Vg;  Vg_crust(~is_crust)  = NaN;
        Vg_mantle = Vg;  Vg_mantle(~is_mantle) = NaN;

        vc = Vg_crust(~isnan(Vg_crust));
        vm = Vg_mantle(~isnan(Vg_mantle));
        has_crust = ~isempty(vc);
        has_mantle = ~isempty(vm);

        %% Color limits per zoom level
        if zmax_plot <= 0.5
            if has_crust, clim_crust = [0.8, 3.2]; end
        elseif zmax_plot <= 5
            if has_crust, clim_crust = [0.8, 3.8]; end
        else
            if has_crust,  clim_crust  = [3.5, 4.1]; end
            if has_mantle, clim_mantle = [4.25, 4.75]; end
        end

        %% Three stacked axes
        fig = figure('Position',[50 50 1300 550],'Color','w','Visible','off');

        ax_box = axes('Position',[0.08 0.13 0.72 0.78]); hold on; box on;
        set(ax_box,'YDir','reverse','FontSize',13,'LineWidth',1.2, ...
            'TickDir','out','XMinorTick','on','YMinorTick','on');
        ax_box.Color = 'none';

        ax_mantle = axes('Position',ax_box.Position); hold on;
        set(ax_mantle,'YDir','reverse','Visible','off','Color','none');

        ax_crust = axes('Position',ax_box.Position); hold on;
        set(ax_crust,'YDir','reverse','Visible','off','Color','none');

        linkaxes([ax_box, ax_mantle, ax_crust]);

        xlims = [dist_arc(1), dist_arc(end)];
        ylims = [0, zmax_plot];

        %% Contourf with Brennan-style settings
        if has_crust
            contourf(ax_crust, dist_sect, zmesh, Vg_crust, v_contours, 'EdgeAlpha', 0.1);
            colormap(ax_crust, flipud(turbo(n_colors)));
            caxis(ax_crust, clim_crust);
        end
        if has_mantle
            contourf(ax_mantle, dist_sect, zmesh, Vg_mantle, v_contours, 'EdgeAlpha', 0.1);
            colormap(ax_mantle, flipud(turbo(n_colors)));
            caxis(ax_mantle, clim_mantle);
        end

        set(ax_crust,  'XLim', xlims, 'YLim', ylims);
        set(ax_mantle, 'XLim', xlims, 'YLim', ylims);
        set(ax_box,    'XLim', xlims, 'YLim', ylims);

        %% Interface lines
        zsed_vis = zsed_sect(:); zsed_vis(zsed_vis > zmax_plot) = NaN;
        zmoh_vis = zmoh_sect(:); zmoh_vis(zmoh_vis > zmax_plot) = NaN;
        plot(ax_box, dist_arc, zsed_vis, '-', 'Color',[0.2 0.6 0.2], 'LineWidth',2.5);
        plot(ax_box, dist_arc, zmoh_vis, '-', 'Color',[0.9 0.0 0.75], 'LineWidth',4);

        %% Station markers
        for is = 1:nsta_line
            plot(ax_box, sta_dist(is), ylims(1), 'v', 'MarkerSize',7, ...
                'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',0.8,'Clipping','off');
        end

        axes(ax_box);
        xlabel('Distance (km)','FontSize',14);
        ylabel('Depth (km)','FontSize',14);
        title(sprintf('Line %s: Vs Cross-Section (0-%s)  %s \\rightarrow %s', ...
            line_id, zoom_tags{iz_plot}, names{1}, names{end}), 'FontSize',14,'FontWeight','bold');

        %% Colorbars
        if has_crust
            cb1 = colorbar(ax_crust,'Location','east');
            cb1.Position = [0.82 0.13 0.02 0.78];
            ylabel(cb1,'Crust Vs (km/s)','FontSize',11);
            set(cb1,'LineWidth',1,'TickDirection','out','FontSize',11);
        end
        if has_mantle
            cb2 = colorbar(ax_mantle,'Location','east');
            cb2.Position = [0.87 0.13 0.02 0.78];
            ylabel(cb2,'Mantle Vs (km/s)','FontSize',11);
            set(cb2,'LineWidth',1,'TickDirection','out','FontSize',11);
        end

        fn = sprintf('xsect_%s_%s', zoom_tags{iz_plot}, line_id);
        exportgraphics(fig, fullfile(linedir,[fn '.png']), 'Resolution', 300);
        exportgraphics(fig, fullfile(linedir,[fn '.pdf']), 'ContentType','image','Resolution',300);
        fprintf('Saved %s\n', fn);
        close(fig);
    end
end

fprintf('\nAll cross-sections saved to %s\n', basedir);