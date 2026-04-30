%% a1_1_compile_results.m
%  Load all station inversion results into a common structure
%  Adapted from Brunsvik et al. for SESAME HBI results

run('a0_parameters_setup.m');

%% Build mdls structure from completed stations
mdls = struct('sta', {}, 'lat', [], 'lon', [], ...
    'model', {}, 'fposterior', {}, 'dir', {});

igdsta = 0;
for is = 1:nsta
    sta_name = summary(is).name;
    sta_path = fullfile(results_dir, sta_name);
    if ~exist(sta_path, 'dir'), continue; end

    files = dir(fullfile(sta_path, 'hbi_corr_*.mat'));
    files = files(~contains({files.name}, 'checkpoint'));
    n_chains = length(files);

    if n_chains < 1, continue; end

    %% Check quality: chain count
    if n_chains < desired_chains
        fprintf('  %s: only %d/%d chains, including anyway\n', sta_name, n_chains, desired_chains);
    end

    %% Load first chain to get prior
    tmp = load(fullfile(files(1).folder, files(1).name));
    r = tmp.results;
    if isfield(r,'prior'), prior = r.prior;
    elseif isfield(tmp,'prior'), prior = tmp.prior;
    else, fprintf('  %s: no prior found, skipping\n', sta_name); continue;
    end

    %% Build model summary from median_vs_summary fields
    model = struct();
    model.med_vs = summary(is).med_vs;
    model.l68 = summary(is).l68;
    model.u68 = summary(is).u68;
    model.l95 = summary(is).l95;
    model.u95 = summary(is).u95;
    model.zinterp = summary(is).zinterp;
    model.zmohav = summary(is).med_zmoh;
    model.zmohsig = (summary(is).u68_zmoh - summary(is).l68_zmoh) / 2;
    model.zsedav = summary(is).med_zsed;
    model.zsedsig = (summary(is).u68_zsed - summary(is).l68_zsed) / 2;
    model.best_logL = summary(is).best_logL;
    model.n_chains = n_chains;
    model.n_samples = summary(is).n_samples;

    %% Extract median Vs at coarse depths
    model.Z = zcoarse;
    model.VSav = interp1(summary(is).zinterp, summary(is).med_vs, zcoarse, 'linear', NaN);

    igdsta = igdsta + 1;
    mdls(1).sta{igdsta,1} = sta_name;
    mdls.lat(igdsta,1) = summary(is).lat;
    mdls.lon(igdsta,1) = summary(is).lon;
    mdls.model{igdsta,1} = model;
    mdls.dir{igdsta,1} = sta_path;

    fprintf('%d. %s: %d chains, logL=%.1f, zmoh=%.1f, zsed=%.2f\n', ...
        igdsta, sta_name, n_chains, model.best_logL, model.zmohav, model.zsedav);
end

fprintf('\nCompiled %d stations\n', igdsta);

%% Save
fresults = fullfile(out_dir, sprintf('compiled_results_%s.mat', STAMP));
save(fresults, 'mdls');
fprintf('Saved to %s\n', fresults);