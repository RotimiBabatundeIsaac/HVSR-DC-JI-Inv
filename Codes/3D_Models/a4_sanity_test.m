%% a2_0_view_results.m
%  Quick map-view of station results: Moho, sediment, Vs at selected depths
%  Adapted from Brunsvik et al. for SESAME HBI results

clc; clear;
run('a0_parameters_setup.m');

mdls = load(fresults).mdls;

%% Reorganize structure and extract map-view parameters
nsta_mdl = length(mdls.sta);
models = struct();
zmoh    = zeros(nsta_mdl, 1);
zmohsig = zeros(nsta_mdl, 1);
zsed    = zeros(nsta_mdl, 1);
zsedsig = zeros(nsta_mdl, 1);
z_vs_check = [0.5, 5, 15, 30, 50, 70];
vs      = zeros(nsta_mdl, length(z_vs_check));

for imd = 1:nsta_mdl
    models(imd).sta   = mdls.sta{imd};
    models(imd).lat   = mdls.lat(imd);
    models(imd).lon   = mdls.lon(imd);
    models(imd).model = mdls.model{imd};
    models(imd).dir   = mdls.dir{imd};

    mdl = mdls.model{imd};
    zmoh(imd)    = mdl.zmohav;
    zmohsig(imd) = mdl.zmohsig;
    zsed(imd)    = mdl.zsedav;
    zsedsig(imd) = mdl.zsedsig;

    for ivs = 1:length(z_vs_check)
        [~, iz_near] = min(abs(mdl.Z - z_vs_check(ivs)));
        vs(imd, ivs) = mdl.VSav(iz_near);
    end
end

%% Map-view parameters (Moho, sediment, errors)
fig1 = figure('Position',[50 50 1400 800],'Color','w');
tiledlayout(2,3,'TileSpacing','compact');

titles_row1 = {'Moho Depth (km)', 'Sed Depth (km)', 'Moho Error (km)'};
vals_row1 = {zmoh, zsed, zmohsig};
titles_row2 = {'Sed Error (km)', 'Best logL', 'N Chains'};
best_logL_vec = zeros(nsta_mdl,1);
n_chains_vec = zeros(nsta_mdl,1);
for imd = 1:nsta_mdl
    best_logL_vec(imd) = mdls.model{imd}.best_logL;
    n_chains_vec(imd) = mdls.model{imd}.n_chains;
end
vals_row2 = {zsedsig, best_logL_vec, n_chains_vec};

for ip = 1:3
    nexttile; hold on;
    m_proj('lambert', 'long', map_lon, 'lat', map_lat);
    m_scatter(mdls.lon, mdls.lat, marker_sz, vals_row1{ip}, 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',0.8);
    if have_states
        for s = 1:length(states)
            m_line(states(s).Lon, states(s).Lat, 'color',[0.5 0.5 0.5], 'linewidth',1);
        end
    end
    colorbar; colormap(gca, turbo);
    m_grid('linestyle','none','linewidth',1,'tickdir','out','fontsize',9);
    title(titles_row1{ip}, 'FontSize',12);
end

for ip = 1:3
    nexttile; hold on;
    m_proj('lambert', 'long', map_lon, 'lat', map_lat);
    m_scatter(mdls.lon, mdls.lat, marker_sz, vals_row2{ip}, 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',0.8);
    if have_states
        for s = 1:length(states)
            m_line(states(s).Lon, states(s).Lat, 'color',[0.5 0.5 0.5], 'linewidth',1);
        end
    end
    colorbar; colormap(gca, turbo);
    m_grid('linestyle','none','linewidth',1,'tickdir','out','fontsize',9);
    title(titles_row2{ip}, 'FontSize',12);
end

sgtitle('Station Parameters Overview', 'FontSize',14, 'FontWeight','bold');
exportgraphics(fig1, fullfile(fig_dir, 'map_view_parameters.png'), 'Resolution', 300);
fprintf('Saved map_view_parameters\n');

%% Pseudo-tomography: Vs at selected depths
fig2 = figure('Position',[50 50 1400 800],'Color','w');
tiledlayout(2,3,'TileSpacing','compact');

for ivs = 1:length(z_vs_check)
    nexttile; hold on;
    m_proj('lambert', 'long', map_lon, 'lat', map_lat);
    m_scatter(mdls.lon, mdls.lat, marker_sz, vs(:,ivs), 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',0.8);
    if have_states
        for s = 1:length(states)
            m_line(states(s).Lon, states(s).Lat, 'color',[0.5 0.5 0.5], 'linewidth',1);
        end
    end
    colorbar; colormap(gca, turbo);
    m_grid('linestyle','none','linewidth',1,'tickdir','out','fontsize',9);
    title(sprintf('Vs at %.1f km', z_vs_check(ivs)), 'FontSize',12);
end

sgtitle('Pseudo-Tomography: Median Vs at Depth', 'FontSize',14, 'FontWeight','bold');
exportgraphics(fig2, fullfile(fig_dir, 'pseudo_tomography.png'), 'Resolution', 300);
fprintf('Saved pseudo_tomography\n');