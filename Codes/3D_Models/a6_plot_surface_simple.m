function [] = a6_plot_surface_simple(llminmax, options)
    arguments
        llminmax
        options.stax = []
        options.stay = []
        options.stav = []
        options.stalon = []
        options.stalat = []
        options.xgrid = []
        options.ygrid = []
        options.vgrid = []
        options.fignum = 1
        options.title = []
        options.sectlon = []
        options.sectlat = []
    end
% Map figure with optional surface contourf and station scatter
% Adapted from Brunsvik et al. for SESAME HBI

lonmin = llminmax(1);
lonmax = llminmax(2);
latmin = llminmax(3);
latmax = llminmax(4);

figure(options.fignum); clf; hold on;
set(gcf, 'color', 'white');
m_proj('lambert', 'long', [lonmin, lonmax], 'lat', [latmin, latmax]);

% Station coordinates
if isempty(options.stax)
    if ~isempty(options.stalon)
        [stax, stay] = m_ll2xy(options.stalon, options.stalat);
    else
        stax = [];
        stay = [];
    end
else
    stax = options.stax;
    stay = options.stay;
end
if isempty(options.stav)
    options.stav = ones(size(stax));
end

% Plot surface (mask points far from stations)
if ~isempty(options.xgrid)
    if ~isempty(stax)
        pt_dist_nan = 0.02;
        pt_dist = zeros(size(options.xgrid));
        for ipt = 1:numel(options.xgrid)
            pt_dist(ipt) = min(sqrt((options.xgrid(ipt) - stax).^2 + ...
                (options.ygrid(ipt) - stay).^2));
        end
        options.vgrid(pt_dist > pt_dist_nan) = nan;
    end
    contourf(options.xgrid, options.ygrid, options.vgrid, 15, 'LineStyle','none');
end

% State boundaries
try
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    for s = 1:length(states)
        m_line(states(s).Lon, states(s).Lat, 'LineWidth', 1, 'color', [0 0 0]);
    end
catch
    try
        [latbord, lonbord] = borders('states');
        for iplace = 1:length(lonbord)
            m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth', 1, 'color', [0 0 0]);
        end
    catch
    end
end

m_grid('box','fancy','linestyle','-','gridcolor', 0.5*[1 1 1], 'backcolor', [0.3 0.75 1]);

if ~isempty(options.title)
    title(options.title, 'fontweight', 'normal');
end

% Plot stations
if ~isempty(stax)
    scatter(stax, stay, 30, options.stav, 'filled', 'MarkerEdgeColor', 'k');
end
cbar = colorbar();
cbar.Label.String = 'Vs';
colormap(turbo);
try
    caxis([min(options.stav), max(options.stav)]);
catch
end

% Cross-section lines
if ~isempty(options.sectlon)
    for isect = 1:size(options.sectlon, 2)
        [sectx, secty] = m_ll2xy(options.sectlon(:,isect), options.sectlat(:,isect));
        plot(sectx, secty, 'k', 'Linewidth', 3);
        section_letter = char(64 + isect);
        text(sectx(1), secty(1), section_letter, 'fontsize', 22, 'color', 'r');
        text(sectx(end), secty(end), [section_letter ''''], 'fontsize', 22, 'color', 'r');
    end
end
end