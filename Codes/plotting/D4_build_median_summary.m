%Build model summary, find the moho and bedrock depth depth 

clear; close all;

results_dir = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/Results_and_Figures/real_data_correlated_results';
datapack_path = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/DATA/sta_datapack.mat';
load(datapack_path);

sta_dirs = dir(results_dir);
sta_dirs = sta_dirs([sta_dirs.isdir] & ~ismember({sta_dirs.name},{'.','..'}));
sta_dirs = sta_dirs(cellfun(@(x) length(x) <= 4 && any(x(1) == 'DEW'), {sta_dirs.name}));

summary = struct();
count = 0;

for is = 1:length(sta_dirs)
    sta_name = sta_dirs(is).name;
    chain_fp = fullfile(results_dir, sta_name, sprintf('hbi_corr_*%s*.mat', sta_name));
    files = dir(chain_fp);
    if isempty(files)
        chain_fp = fullfile(results_dir, sta_name, 'hbi_corr_*.mat');
        files = dir(chain_fp);
    end
    if isempty(files), continue; end

    n_chains = length(files);
    all_theta = []; all_logL = [];
    prior = [];
    for k = 1:n_chains
        tmp = load(fullfile(files(k).folder, files(k).name));
        r = tmp.results;
        all_theta = [all_theta; r.samples_theta];
        all_logL  = [all_logL; r.samples_logL(:)];
        if k == 1
            if isfield(r,'prior'), prior = r.prior;
            elseif isfield(tmp,'prior'), prior = tmp.prior; end
        end
    end
    if isempty(prior), continue; end

    idx_h = prior.idx.h; idx_vs = prior.idx.vs;
    h_all = all_theta(:,idx_h); vs_all = all_theta(:,idx_vs);
    zmax = sum(prior.bounds.upper(idx_h));
    vsmax = max(prior.bounds.upper(idx_vs));
    zmax_actual = max(sum(h_all,2));

    dz = 1/1000;
    zinterp = 0:dz:zmax;
    nz = length(zinterp);
    Nmods = size(all_theta,1);
    L_norm = exp(all_logL - max(all_logL));

    vs_mod = nan(Nmods, nz);
    for imod = 1:Nmods
        hmod = [h_all(imod,:)'; zmax - sum(h_all(imod,:),2)];
        vsmod = vs_all(imod,:)';
        Igood = all_theta(imod,:) >= prior.bounds.lower' & all_theta(imod,:) <= prior.bounds.upper';
        Igood_vs = Igood(idx_vs); Igood_h = [Igood(idx_h), true];
        prior_mod = Igood_h(:) .* Igood_vs(:);
        modn = [hmod, vsmod*0, vsmod, prior_mod];
        mod_int = layerizemod_interp(modn, zinterp);
        vs_mod(imod,:) = mod_int.vs;
    end

    med_vs = nan(1,nz); l68 = nan(1,nz); u68 = nan(1,nz);
    l95 = nan(1,nz); u95 = nan(1,nz);
    for iz = 1:nz
        vord = sort(vs_mod(:,iz));
        med_vs(iz) = vord(round(0.5*Nmods));
        l68(iz) = vord(max(1, round(0.16*Nmods)));
        u68(iz) = vord(min(Nmods, round(0.84*Nmods)));
        l95(iz) = vord(max(1, round(0.025*Nmods)));
        u95(iz) = vord(min(Nmods, round(0.975*Nmods)));
    end

    z_interfaces = cumsum(h_all, 2);

    zmoh_samples = nan(Nmods, 1);
    zsed_samples = nan(Nmods, 1);
    for imod = 1:Nmods
        vs_vals = vs_all(imod,:);
        dv = diff(vs_vals);

        %moho_candidates = find(vs_vals(1:end-1) < 4.0 & vs_vals(2:end) >= 3.0 & dv > 0);
        moho_candidates = find(vs_vals(1:end-1) >= 3.0 & vs_vals(2:end) >= 4.0 & dv > 0);
        if ~isempty(moho_candidates)
            [~, ibest_moh] = max(dv(moho_candidates));
            imoh = moho_candidates(ibest_moh);
            zmoh_samples(imod) = z_interfaces(imod, imoh);
        end

        % shallow = min(5, length(dv));
        % [~, ised] = max(dv(1:shallow));
        % zsed_samples(imod) = z_interfaces(imod, ised);

            shallow = min(5, length(dv));
            sed_candidates = find(vs_vals(1:shallow) < 2.0 & vs_vals(2:shallow+1) >= 2.0 & dv(1:shallow) > 0);
            if ~isempty(sed_candidates)
                [~, ibest_sed] = max(dv(sed_candidates));
                ised = sed_candidates(ibest_sed);
                zsed_samples(imod) = z_interfaces(imod, ised);
            else
                zsed_samples(imod) = NaN;
            end


    end

    valid_moh = ~isnan(zmoh_samples);
    valid_sed = ~isnan(zsed_samples);
    med_zmoh = median(zmoh_samples(valid_moh));
    l68_zmoh = prctile(zmoh_samples(valid_moh), 16);
    u68_zmoh = prctile(zmoh_samples(valid_moh), 84);
    med_zsed = median(zsed_samples(valid_sed));
    l68_zsed = prctile(zsed_samples(valid_sed), 16);
    u68_zsed = prctile(zsed_samples(valid_sed), 84);

    if sum(valid_moh) < 0.5*Nmods
        fprintf('  Warning: %s Moho detected in only %d/%d samples\n', sta_name, sum(valid_moh), Nmods);
    end

    [best_logL, ibest] = max(all_logL);

    sta_idx = find(strcmp(cellfun(@(x) strrep(x,'Z9.',''), {sta.name}, 'UniformOutput', false), sta_name));
    if isempty(sta_idx), continue; end
    sta_lat = sta(sta_idx(1)).lat;
    sta_lon = sta(sta_idx(1)).lon;

    count = count + 1;
    summary(count).name = sta_name;
    summary(count).lat = sta_lat;
    summary(count).lon = sta_lon;
    summary(count).zinterp = zinterp;
    summary(count).med_vs = med_vs;
    summary(count).l68 = l68;
    summary(count).u68 = u68;
    summary(count).l95 = l95;
    summary(count).u95 = u95;
    summary(count).zmax_actual = zmax_actual;
    summary(count).vsmax = vsmax;
    summary(count).best_logL = best_logL;
    summary(count).n_chains = n_chains;
    summary(count).n_samples = Nmods;
    summary(count).med_zsed = med_zsed;
    summary(count).l68_zsed = l68_zsed;
    summary(count).u68_zsed = u68_zsed;
    summary(count).med_zmoh = med_zmoh;
    summary(count).l68_zmoh = l68_zmoh;
    summary(count).u68_zmoh = u68_zmoh;
    fprintf('%d. %s: %d chains, %d samples, best logL=%.1f, zsed=%.2f km, zmoh=%.1f km\n', ...
        count, sta_name, n_chains, Nmods, best_logL, med_zsed, med_zmoh);
end

outfile = fullfile(results_dir, 'median_vs_summary.mat');
save(outfile, 'summary', '-v7.3');
fprintf('\nSaved %d stations to %s\n', count, outfile);