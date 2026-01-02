%==========================================================================
% build_station_datapack.m  (DEBUGGED v2 - with diagnostics)
%==========================================================================
clear; clc;

%% ---------------------- USER PATHS -----------------------
ROOT = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/DATA';

grv_mat     = fullfile(ROOT, 'Gr_Vel', 'raytomo_grv_3_25s_ZZ.mat');
phv_mat     = fullfile(ROOT, 'Ph_Vel', 'raytomo_9_35s_1wl_ZZ_0S_LRT.mat');
hvstats_dir = fullfile(ROOT, 'HV_ellip_stats');
hvsr_dir    = fullfile(ROOT, 'HVSR');

Npixel_dist = 2;

sig_model.phase_rel = 0.02;  % 2% relative
sig_model.group_rel = 0.03;  % 3% relative
sig_model.min_abs   = 0.02;  % km/s floor

out_mat = fullfile(ROOT, 'sta_datapack.mat');

DEBUG = true;  % Set to false to suppress verbose output

%% ---------------------- CHECKS ---------------------------
assert(exist(grv_mat,'file')==2, 'Missing GRV mat: %s', grv_mat);
assert(exist(phv_mat,'file')==2, 'Missing PHV mat: %s', phv_mat);
assert(exist(hvstats_dir,'dir')==7, 'Missing HV stats dir: %s', hvstats_dir);
assert(exist(hvsr_dir,'dir')==7, 'Missing HVSR dir: %s', hvsr_dir);

%% ---------------------- LOAD GRV -------------------------
fprintf('Loading GRV: %s\n', grv_mat);
G = load(grv_mat);
para_grv = G.parameters;
raytomo_grv = G.raytomo;
gridsize_grv = para_grv.gridsize;

if DEBUG
    fprintf('  GRV file contains fields: %s\n', strjoin(fieldnames(G), ', '));
    fprintf('  raytomo_grv has %d period entries\n', numel(raytomo_grv));
    fprintf('  gridsize_grv = %.4f deg\n', gridsize_grv);
end

% --- GRV grid extraction (inline) ---
[latgrid_grv, longrid_grv, grv_grid_source] = extract_latlon_grids(G, 'GRV');
if DEBUG
    fprintf('  GRV grid source: %s\n', grv_grid_source);
    fprintf('  GRV lat range: [%.2f, %.2f], lon range: [%.2f, %.2f]\n', ...
        min(latgrid_grv(:)), max(latgrid_grv(:)), ...
        min(longrid_grv(:)), max(longrid_grv(:)));
end

periods_grv = [raytomo_grv.period];
nper_grv = numel(periods_grv);

station_list_file = para_grv.station_list;
assert(exist(station_list_file,'file')==2, 'Station list not found: %s', station_list_file);

[net, name, lat, lon, dep] = textread(station_list_file, '%s %s %f %f %f'); %#ok<DTXTRD>
nsta = numel(name);
sta_full = cell(nsta,1);
for i = 1:nsta
    sta_full{i} = [net{i}, '.', name{i}];
end

if DEBUG
    fprintf('  Loaded %d stations from: %s\n', nsta, station_list_file);
    fprintf('  Station lat range: [%.2f, %.2f], lon range: [%.2f, %.2f]\n', ...
        min(lat), max(lat), min(lon), max(lon));
end

%% ---------------------- LOAD PHV -------------------------
fprintf('\nLoading PHV: %s\n', phv_mat);
P = load(phv_mat);
raytomo_phv = P.raytomo;
periods_phv = [raytomo_phv.period];
nper_phv = numel(periods_phv);

if DEBUG
    fprintf('  PHV file contains fields: %s\n', strjoin(fieldnames(P), ', '));
    fprintf('  raytomo_phv has %d period entries\n', numel(raytomo_phv));
end

% --- PHV grid extraction ---
[latgrid_phv, longrid_phv, phv_grid_source] = extract_latlon_grids(P, 'PHV');
if DEBUG
    fprintf('  PHV grid source: %s\n', phv_grid_source);
    fprintf('  PHV lat range: [%.2f, %.2f], lon range: [%.2f, %.2f]\n', ...
        min(latgrid_phv(:)), max(latgrid_phv(:)), ...
        min(longrid_phv(:)), max(longrid_phv(:)));
end

phv_field = detect_phv_field(raytomo_phv);
if isempty(phv_field)
    warning('No obvious PHV field found. Phase curves will be NaN.');
else
    fprintf('  Detected PHV field: %s\n', phv_field);
end

gridsize_phv = gridsize_grv;
if isfield(P,'parameters') && isfield(P.parameters,'gridsize') && ~isempty(P.parameters.gridsize)
    gridsize_phv = P.parameters.gridsize;
end

%% ---------------------- BUILD STA STRUCT ------------------
sta = repmat(struct(), nsta, 1);

GVmap0 = raytomo_grv(1).GV;
phvMap0 = [];
if ~isempty(phv_field) && isfield(raytomo_phv(1), phv_field)
    phvMap0 = raytomo_phv(1).(phv_field);
end

if DEBUG
    fprintf('\n  GV map size: %s\n', mat2str(size(GVmap0)));
    if ~isempty(phvMap0)
        fprintf('  PHV map size: %s\n', mat2str(size(phvMap0)));
    end
end

fprintf('\nSampling GRV and PHV at station locations...\n');
grv_sample_fails = {};
phv_sample_fails = {};

for is = 1:nsta
    sta(is).name = sta_full{is};
    sta(is).lat  = lat(is);
    sta(is).lon  = lon(is);
    sta(is).dep  = dep(is);

    % Distance grid for GRV
    [rkm_grv, grv_ok] = build_distance_grid(lat(is), lon(is), latgrid_grv, longrid_grv, size(GVmap0), DEBUG && is==1);
    if ~grv_ok
        grv_sample_fails{end+1} = sta_full{is}; %#ok<AGROW>
    end

    % Distance grid for PHV
    rkm_phv = [];
    phv_ok = true;
    if ~isempty(phvMap0)
        [rkm_phv, phv_ok] = build_distance_grid(lat(is), lon(is), latgrid_phv, longrid_phv, size(phvMap0), false);
        if ~phv_ok
            phv_sample_fails{end+1} = sta_full{is}; %#ok<AGROW>
        end
    end

    % Group velocity
    gv_vec = nan(nper_grv, 1);
    if grv_ok
        for ip = 1:nper_grv
            gv_vec(ip) = sample_nearest_good_pixel(raytomo_grv(ip).GV, rkm_grv, gridsize_grv, Npixel_dist);
        end
    end
    sta(is).ugr.T   = periods_grv(:);
    sta(is).ugr.obs = gv_vec;
    sta(is).ugr.sig = build_disp_sig(gv_vec, sig_model.group_rel, sig_model.min_abs);

    % Phase velocity - handle both 2D grid and 1D pixel-list formats
    phv_vec = nan(nper_phv, 1);
    if ~isempty(phv_field)
        for ip = 1:nper_phv
            if ~isfield(raytomo_phv(ip), phv_field)
                continue;
            end
            phv_map = raytomo_phv(ip).(phv_field);
            
            % Check if PHV is stored as 1D pixel list
            if isvector(phv_map)
                % Pixel-list format: phv(i) at location (latgrid_phv(i), longrid_phv(i))
                phv_vec(ip) = sample_from_pixel_list(lat(is), lon(is), ...
                    latgrid_phv, longrid_phv, phv_map, gridsize_phv, Npixel_dist);
            else
                % Standard 2D grid format
                if phv_ok && ~isempty(rkm_phv)
                    phv_vec(ip) = sample_nearest_good_pixel(phv_map, rkm_phv, gridsize_phv, Npixel_dist);
                end
            end
        end
    end
    sta(is).cph.T   = periods_phv(:);
    sta(is).cph.obs = phv_vec;
    sta(is).cph.sig = build_disp_sig(phv_vec, sig_model.phase_rel, sig_model.min_abs);
end

if ~isempty(grv_sample_fails)
    fprintf('  WARNING: GRV grid mismatch for %d stations\n', numel(grv_sample_fails));
end
if ~isempty(phv_sample_fails)
    fprintf('  WARNING: PHV grid mismatch for %d stations\n', numel(phv_sample_fails));
end

%% ---------------------- ELLIP (HV_STATS) -----------------
fprintf('\nLoading ellipticity hvstats from: %s\n', hvstats_dir);
ellip_found = {};
ellip_missing = {};

for is = 1:nsta
    [T, obs, sig, found_file] = load_hvstats_for_station(hvstats_dir, sta(is).name);
    sta(is).ellip.T   = T;
    sta(is).ellip.obs = obs;
    sta(is).ellip.sig = sig;
    if ~isempty(found_file)
        ellip_found{end+1} = sta(is).name; %#ok<AGROW>
    else
        ellip_missing{end+1} = sta(is).name; %#ok<AGROW>
    end
end

if DEBUG && ~isempty(ellip_missing)
    fprintf('  Ellip files missing for %d stations (first 5): %s\n', ...
        numel(ellip_missing), strjoin(ellip_missing(1:min(5,end)), ', '));
end

%% ---------------------- HVSR ------------------------------
fprintf('\nLoading HVSR txt from: %s\n', hvsr_dir);
hvsr_found = {};
hvsr_missing = {};

for is = 1:nsta
    [f, obs, sig, found_file] = load_hvsr_for_station(hvsr_dir, sta(is).name, DEBUG && is==1);
    sta(is).hvsr.f   = f;
    sta(is).hvsr.obs = obs;
    sta(is).hvsr.sig = sig;
    if ~isempty(found_file)
        hvsr_found{end+1} = sta(is).name; %#ok<AGROW>
    else
        hvsr_missing{end+1} = sta(is).name; %#ok<AGROW>
    end
end

if DEBUG && ~isempty(hvsr_missing)
    fprintf('  HVSR files missing for %d stations (first 5): %s\n', ...
        numel(hvsr_missing), strjoin(hvsr_missing(1:min(5,end)), ', '));
end

%% ---------------------- SAVE ------------------------------
save(out_mat, 'sta', 'grv_mat', 'phv_mat', 'hvstats_dir', 'hvsr_dir', 'station_list_file');
fprintf('\nSaved station datapack: %s\n', out_mat);

%% ---------------------- REPORT ----------------------------
n_hvsr  = sum(arrayfun(@(s) ~isempty(s.hvsr.f), sta));
n_ellip = sum(arrayfun(@(s) ~isempty(s.ellip.T), sta));
n_grv   = sum(arrayfun(@(s) any(isfinite(s.ugr.obs)), sta));
n_phv   = sum(arrayfun(@(s) any(isfinite(s.cph.obs)), sta));

fprintf('\n========== SANITY REPORT (out of %d stations) ==========\n', nsta);
fprintf('  HVSR loaded for:   %3d  (%.1f%%)\n', n_hvsr, 100*n_hvsr/nsta);
fprintf('  Ellip loaded for:  %3d  (%.1f%%)\n', n_ellip, 100*n_ellip/nsta);
fprintf('  GRV sampled for:   %3d  (%.1f%%)\n', n_grv, 100*n_grv/nsta);
fprintf('  PHV sampled for:   %3d  (%.1f%%)\n', n_phv, 100*n_phv/nsta);
fprintf('==========================================================\n');

% Quick data quality check
if n_grv < nsta*0.5
    warning('Less than 50%% of stations have valid GRV data!');
end
if n_phv < nsta*0.5
    warning('Less than 50%% of stations have valid PHV data!');
end

%% ========================= HELPERS =========================

function [latgrid, longrid, source] = extract_latlon_grids(S, label)
    % Try multiple possible field names for grid coordinates
    source = 'unknown';
    latgrid = []; longrid = [];
    
    % Priority 1: xi/yi (common in raytomo output)
    if isfield(S,'xi') && isfield(S,'yi')
        [latgrid, longrid] = infer_latlon_grids(S.xi, S.yi);
        source = 'xi/yi';
        return;
    end
    
    % Priority 2: xnode/ynode
    if isfield(S,'xnode') && isfield(S,'ynode')
        [latgrid, longrid] = infer_latlon_grids(S.xnode, S.ynode);
        source = 'xnode/ynode';
        return;
    end
    
    % Priority 3: lalim/lolim with gridsize (reconstruct grid)
    if isfield(S,'parameters')
        p = S.parameters;
        if isfield(p,'lalim') && isfield(p,'lolim') && isfield(p,'gridsize')
            latgrid = p.lalim(1):p.gridsize:p.lalim(2);
            longrid = p.lolim(1):p.gridsize:p.lolim(2);
            source = 'lalim/lolim reconstruction';
            return;
        end
    end
    
    error('%s mat missing grid fields. Available fields: %s', label, strjoin(fieldnames(S), ', '));
end

function [latgrid, longrid] = infer_latlon_grids(A, B)
    a = A(:); b = B(:);
    a = a(isfinite(a)); b = b(isfinite(b));
    
    if isempty(a) || isempty(b)
        warning('Empty grids; defaulting to B=lat, A=lon.');
        latgrid = B; longrid = A;
        return;
    end
    
    amin = min(a); amax = max(a);
    bmin = min(b); bmax = max(b);

    % Check for clearly negative values (Western Hemisphere longitude)
    % Longitude in W. Hemisphere: typically -180 to -20
    % Latitude: typically -90 to +90, but SE US is +25 to +50
    
    a_is_negative_lon = (amax < -20);  % All values negative and < -20 => W. Hemisphere lon
    b_is_negative_lon = (bmax < -20);
    
    % If one is clearly negative longitude, use that
    if b_is_negative_lon && ~a_is_negative_lon
        longrid = B; latgrid = A;
        return;
    elseif a_is_negative_lon && ~b_is_negative_lon
        longrid = A; latgrid = B;
        return;
    end
    
    % Check for clearly positive latitude (Northern Hemisphere, non-polar)
    % SE US latitudes: ~25 to 50
    a_is_positive_lat = (amin > 0) && (amax < 90) && (amax < 60);
    b_is_positive_lat = (bmin > 0) && (bmax < 90) && (bmax < 60);
    
    if a_is_positive_lat && ~b_is_positive_lat
        latgrid = A; longrid = B;
        return;
    elseif b_is_positive_lat && ~a_is_positive_lat
        latgrid = B; longrid = A;
        return;
    end
    
    % Fallback: use absolute value ranges
    % Longitude typically has larger absolute values or range
    a_range = amax - amin;
    b_range = bmax - bmin;
    
    if abs(mean(b)) > abs(mean(a))
        longrid = B; latgrid = A;
    elseif abs(mean(a)) > abs(mean(b))
        longrid = A; latgrid = B;
    elseif a_range > b_range
        longrid = A; latgrid = B;
    else
        longrid = B; latgrid = A;
    end
end

function [rkm, ok] = build_distance_grid(sta_lat, sta_lon, latgrid, longrid, mapSize, verbose)
    ok = true;
    
    if isvector(latgrid) && isvector(longrid)
        latv = latgrid(:);
        lonv = longrid(:);
        [LON, LAT] = meshgrid(lonv, latv);
    else
        LAT = latgrid;
        LON = longrid;
    end

    % Check size compatibility
    if ~isequal(size(LAT), mapSize)
        if isequal(size(LAT'), mapSize)
            LAT = LAT'; LON = LON';
        elseif isequal([numel(latgrid), numel(longrid)], mapSize)
            % Grids are vectors but meshgrid gave wrong orientation
            [LON, LAT] = meshgrid(longrid(:), latgrid(:));
        elseif isequal([numel(longrid), numel(latgrid)], mapSize)
            [LAT, LON] = meshgrid(latgrid(:), longrid(:));
        else
            if verbose
                fprintf('    WARNING: Grid size %s does not match map size %s\n', ...
                    mat2str(size(LAT)), mat2str(mapSize));
            end
            ok = false;
            rkm = nan(mapSize);
            return;
        end
    end

    rkm = haversine_km(sta_lat, sta_lon, LAT, LON);
    
    if verbose
        fprintf('    Distance grid built: size=%s, min_dist=%.2f km\n', ...
            mat2str(size(rkm)), min(rkm(:)));
    end
end

function dkm = haversine_km(lat1, lon1, lat2, lon2)
    R = 6371.0;
    phi1 = deg2rad(lat1);
    phi2 = deg2rad(lat2);
    dphi = deg2rad(lat2 - lat1);
    dlmb = deg2rad(lon2 - lon1);

    a = sin(dphi/2).^2 + cos(phi1).*cos(phi2).*sin(dlmb/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    dkm = R * c;
end

function val = sample_nearest_good_pixel(Vmap, rkm, gridsize_deg, Npixel_dist)
    % Handle size mismatch
    if ~isequal(size(rkm), size(Vmap))
        if isequal(size(rkm'), size(Vmap))
            rkm = rkm';
        else
            val = NaN;
            return;
        end
    end

    rkm_nan = rkm;
    rkm_nan(isnan(Vmap) | Vmap <= 0) = NaN;  % Also reject non-physical velocities

    [rmin, imin] = min(rkm_nan(:));
    if isempty(rmin) || isnan(rmin)
        val = NaN;
        return;
    end

    pix_km = 111.195 * gridsize_deg;
    if rmin / pix_km > Npixel_dist
        val = NaN;
        return;
    end

    val = Vmap(imin);
end

function val = sample_from_pixel_list(sta_lat, sta_lon, lat_coords, lon_coords, values, gridsize_deg, Npixel_dist)
    % Sample from a 1D pixel-list format where values(i) is at (lat_coords(i), lon_coords(i))
    % Flatten all inputs
    lat_vec = lat_coords(:);
    lon_vec = lon_coords(:);
    val_vec = values(:);
    
    % Compute distances from station to all pixels
    rkm = haversine_km(sta_lat, sta_lon, lat_vec, lon_vec);
    
    % Mask invalid values
    rkm(isnan(val_vec) | val_vec <= 0) = NaN;
    
    % Find nearest valid pixel
    [rmin, imin] = min(rkm);
    if isempty(rmin) || isnan(rmin)
        val = NaN;
        return;
    end
    
    % Check distance threshold
    pix_km = 111.195 * gridsize_deg;
    if rmin / pix_km > Npixel_dist
        val = NaN;
        return;
    end
    
    val = val_vec(imin);
end

function phv_field = detect_phv_field(raytomo_phv)
    % Detect the phase velocity field, preferring 2D gridded fields over 1D vectors
    phv_field = '';
    
    % Priority 1: 2D gridded fields (preferred)
    candidates_2d = {'GV', 'phv_iso_2d', 'c_iso', 'Ciso', 'cRayl', 'c_rayl'};
    for k = 1:numel(candidates_2d)
        if isfield(raytomo_phv, candidates_2d{k})
            testval = raytomo_phv(1).(candidates_2d{k});
            % Must be 2D (not a vector) and contain valid data
            if ~isempty(testval) && ~isvector(testval) && any(isfinite(testval(:)))
                phv_field = candidates_2d{k};
                return;
            end
        end
    end
    
    % Priority 2: Any field with valid data (fallback)
    candidates_any = {'phv', 'c', 'C', 'PHV'};
    for k = 1:numel(candidates_any)
        if isfield(raytomo_phv, candidates_any{k})
            testval = raytomo_phv(1).(candidates_any{k});
            if ~isempty(testval) && any(isfinite(testval(:)))
                phv_field = candidates_any{k};
                return;
            end
        end
    end
end

function sig = build_disp_sig(v, rel, min_abs)
    sig = rel .* abs(v);
    sig(sig < min_abs) = min_abs;
    sig(~isfinite(v)) = NaN;
end

function [T, obs, sig, found_file] = load_hvstats_for_station(hvstats_dir, sta_full)
    T = []; obs = []; sig = []; found_file = '';
    parts = strsplit(sta_full, '.'); suf = parts{end};

    cand = {
        fullfile(hvstats_dir, [sta_full, '_hvstats.mat'])
        fullfile(hvstats_dir, [suf, '_hvstats.mat'])
    };

    for i = 1:numel(cand)
        if exist(cand{i}, 'file') == 2
            found_file = cand{i};
            break;
        end
    end

    % Fuzzy match if exact match fails
    if isempty(found_file)
        files = dir(fullfile(hvstats_dir, '*_hvstats.mat'));
        key_full = normalize_name(sta_full);
        key_suf  = normalize_name(suf);
        for i = 1:numel(files)
            base = strrep(files(i).name, '_hvstats.mat', '');
            bkey = normalize_name(base);
            if strcmp(bkey, key_full) || strcmp(bkey, key_suf)
                found_file = fullfile(hvstats_dir, files(i).name);
                break;
            end
        end
    end

    if isempty(found_file), return; end

    S = load(found_file);
    if ~isfield(S,'hv_stats'), return; end
    hv_stats = S.hv_stats;

    if ~isfield(hv_stats,'periods'), return; end
    T = hv_stats.periods(:);

    if isfield(hv_stats,'hv_med')
        obs = hv_stats.hv_med(:);
    elseif isfield(hv_stats,'hv_mean')
        obs = hv_stats.hv_mean(:);
    else
        T = []; return;
    end

    if isfield(hv_stats,'hv_std')
        sig = hv_stats.hv_std(:);
    else
        sig = 0.1 * abs(obs);  % Default 10% uncertainty
    end
end

function [f, obs, sig, found_file] = load_hvsr_for_station(hvsr_dir, sta_full, verbose)
    f = []; obs = []; sig = []; found_file = '';
    parts = strsplit(sta_full, '.'); sta = parts{end};

    cand = {
        fullfile(hvsr_dir, [sta, '_meanstd.txt'])
        fullfile(hvsr_dir, [sta_full, '_meanstd.txt'])
        fullfile(hvsr_dir, [sta, '.txt'])
        fullfile(hvsr_dir, [sta_full, '.txt'])
        fullfile(hvsr_dir, [sta, '_hvsr.txt'])
    };

    for i = 1:numel(cand)
        if exist(cand{i}, 'file') == 2
            found_file = cand{i};
            break;
        end
    end
    if isempty(found_file), return; end

    % Read entire file, skip header lines starting with #
    % Format: "0.20, 0.94, 0.11" (comma-separated)
    fid = fopen(found_file,'r');
    if fid < 0
        warning('Could not open HVSR file: %s', found_file);
        return;
    end
    
    raw = textscan(fid, '%f %f %f', 'Delimiter', ',', 'CommentStyle', '#');
    fclose(fid);
    
    if isempty(raw{1})
        if verbose
            fprintf('    HVSR file %s: no valid data rows\n', found_file);
        end
        found_file = '';
        return;
    end
    
    f   = raw{1};
    obs = raw{2};
    sig = raw{3};
    
    % Remove rows with NaN
    valid = isfinite(f) & isfinite(obs) & isfinite(sig);
    f = f(valid); obs = obs(valid); sig = sig(valid);
    
    if isempty(f)
        found_file = '';
        return;
    end

    % Handle duplicate frequencies by averaging
    [fu, ~, ic] = unique(f);
    if numel(fu) < numel(f)
        obs = accumarray(ic, obs, [], @mean);
        sig = accumarray(ic, sig, [], @mean);
        f = fu;
    end
    
    if verbose
        fprintf('    HVSR loaded: %d points, f=[%.2f - %.2f] Hz\n', numel(f), min(f), max(f));
    end
end

function key = normalize_name(s)
    key = upper(regexprep(strtrim(s), '[\.\_\-\s]', ''));
end