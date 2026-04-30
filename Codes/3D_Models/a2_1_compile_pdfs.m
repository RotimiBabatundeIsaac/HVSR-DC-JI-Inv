%% a2_1_compile_pdfs.m
%  Compile posterior PDFs from all station results for surface inversion


run('a0_parameters_setup.m');
addpath('/Volumes/Drive/AS-Filer/EES/jbrussel/SharedData/SESAME/calc_HBI_posterior_jbr');
%% Parameters
nkernel = 100;
widthkernel_vel = 0.03;
widthkernel_parm = 0.02;

%% Output file
fpdfs = fullfile(out_dir, sprintf('compiled_pdfs_%s.mat', STAMP));

%% Initiate PDF structures
pdfs_empty = struct('sta', {}, 'mm', {}, 'pm', {});

pdfs_allparm = struct();
indiv_parameters = ["zsed", "zmoh"];
for iparam = 1:length(indiv_parameters)
    fn = indiv_parameters(iparam);
    pdfs_allparm(1).(fn) = pdfs_empty;
end

pdfs_allparm.vs = {};
for iz = 1:nz
    pdfs_allparm(1).vs{iz} = pdfs_empty;
end

%% Loop over stations and build PDFs
for is = 1:nsta
    sta_name = summary(is).name;

    %% Load chain files for this station
    sta_path = fullfile(results_dir, sta_name);
    files = dir(fullfile(sta_path, 'hbi_corr_*.mat'));
    files = files(~contains({files.name}, 'checkpoint'));
    if isempty(files), continue; end

    all_theta = []; all_logL = [];
    prior = [];
    for k = 1:length(files)
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

    idx_h = prior.idx.h;
    idx_vs = prior.idx.vs;
    h_all = all_theta(:, idx_h);
    vs_all = all_theta(:, idx_vs);
    Nmods = size(all_theta, 1);

    %% Moho PDF from posterior samples
    z_interfaces = cumsum(h_all, 2);
    zmoh_samples = nan(Nmods, 1);
    zsed_samples = nan(Nmods, 1);

    for imod = 1:Nmods
        vs_vals = vs_all(imod,:);
        dv = diff(vs_vals);

        moho_candidates = find(vs_vals(1:end-1) >= 3.0 & vs_vals(2:end) >= 4.0 & dv > 0);
        if ~isempty(moho_candidates)
            [~, ibest_moh] = max(dv(moho_candidates));
            imoh = moho_candidates(ibest_moh);
            zmoh_samples(imod) = z_interfaces(imod, imoh);
        end

        shallow = min(5, length(dv));
        sed_candidates = find(vs_vals(1:shallow) < 2.0 & vs_vals(2:shallow+1) >= 2.0 & dv(1:shallow) > 0);
        if ~isempty(sed_candidates)
            [~, ibest_sed] = max(dv(sed_candidates));
            ised = sed_candidates(ibest_sed);
            zsed_samples(imod) = z_interfaces(imod, ised);
        end
    end

    %% KDE for Moho depth
    valid_moh = zmoh_samples(~isnan(zmoh_samples));
    if length(valid_moh) > 10
        dsamp_max = max(valid_moh) - min(valid_moh);
        dsamp_max = max(dsamp_max, 0.0001);
        wk = widthkernel_parm * dsamp_max;
        [pdfm, mm] = ksdensity(valid_moh, 'width', wk, 'NumPoints', nkernel);
        pdfs_allparm(is).zmoh(1).mm = mm';
        pdfs_allparm(is).zmoh(1).pm = pdfm';
        pdfs_allparm(is).zmoh(1).sta = sta_name;
    end

    %% KDE for sediment depth
    valid_sed = zsed_samples(~isnan(zsed_samples));
    if length(valid_sed) > 10
        dsamp_max = max(valid_sed) - min(valid_sed);
        dsamp_max = max(dsamp_max, 0.0001);
        wk = widthkernel_parm * dsamp_max;
        [pdfm, mm] = ksdensity(valid_sed, 'width', wk, 'NumPoints', nkernel);
        pdfs_allparm(is).zsed(1).mm = mm';
        pdfs_allparm(is).zsed(1).pm = pdfm';
        pdfs_allparm(is).zsed(1).sta = sta_name;
    end

    %% Interpolate Vs profiles to coarse depth grid
    zmax = sum(prior.bounds.upper(idx_h));
    dz_fine = 1/1000;
    zinterp_fine = 0:dz_fine:zmax;
    nz_fine = length(zinterp_fine);

    vs_mod = nan(Nmods, nz_fine);
    for imod = 1:Nmods
        hmod = [h_all(imod,:)'; zmax - sum(h_all(imod,:),2)];
        vsmod = vs_all(imod,:)';
        Igood = all_theta(imod,:) >= prior.bounds.lower' & all_theta(imod,:) <= prior.bounds.upper';
        Igood_vs = Igood(idx_vs);
        Igood_h = [Igood(idx_h), true];
        prior_mod = Igood_h(:) .* Igood_vs(:);
        modn = [hmod, vsmod*0, vsmod, prior_mod];
        mod_int = layerizemod_interp(modn, zinterp_fine);
        vs_mod(imod,:) = mod_int.vs;
    end

    %% KDE for velocity at each coarse depth
    dv_max = max(vs_mod(:)) - min(vs_mod(:));
    kernel_points = linspace(min(vs_mod(:)), max(vs_mod(:)), nkernel)';

    pdfs_allparm(is).zcoarse = zcoarse;
    for iz = 1:nz
        [~, iz_near] = min(abs(zinterp_fine - zcoarse(iz)));
        samps = vs_mod(:, iz_near);
        samps = samps(~isnan(samps));
        if length(samps) < 10, continue; end

        [pdfm, mm] = ksdensity(samps, kernel_points, 'width', widthkernel_vel);
        pdfs_allparm(is).vs{iz}(1).mm = mm';
        pdfs_allparm(is).vs{iz}(1).pm = pdfm';
        pdfs_allparm(is).vs{iz}(1).sta = sta_name;
    end

    fprintf('%d/%d (%.0f%%) %s: %d chains, %d samples\n', ...
        is, nsta, is/nsta*100, sta_name, length(files), Nmods);
end

save(fpdfs, 'pdfs_allparm', '-v7.3');
fprintf('\nSaved PDFs to %s\n', fpdfs);