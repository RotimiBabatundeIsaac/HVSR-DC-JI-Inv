%% a3_0_surface_inversion.m
%  3-D surface inversion using fminunc with PDF likelihood + roughness penalty
%  Adapted from Brunsvik et al. for SESAME HBI results

clc; clear;
run('a0_parameters_setup.m');

fpdfs = fullfile(out_dir, sprintf('compiled_pdfs_%s.mat', STAMP));
mdls = load(fresults).mdls;

%% Parameters
rough_scale_base = 1e-9;
re_run = true;
max_inv_iterations = 30;
version_surf = 1;

rough_scale_params = struct('zmoh', 0.0015*rough_scale_base, ...
    'zsed', 0.1*rough_scale_base);

% disconts = {"zsed", "zmoh"};
disconts = {"zsed", "zmoh"};
to_invert = disconts;

%% Build inversion list: discontinuities first, then depth slices
to_invert = disconts;
% for inum = 1:nz
%     to_invert{end+1} = inum;
% end

%% Load PDFs
pdf_file = load(fpdfs);
pdfs_all = pdf_file.pdfs_allparm;

nsta_mdl = length(mdls.sta);

%% Set up projected coordinate grid (ndgrid ordering, matching Brennan)
nodes_per_degree = 4;
nx = ceil((llminmax(2) - llminmax(1)) * nodes_per_degree);
ny = ceil((llminmax(4) - llminmax(3)) * nodes_per_degree);
edge_space = 0.2;

DX = max(stax) - min(stax);
DY = max(stay) - min(stay);
xline = linspace(min(stax) - edge_space*DX, max(stax) + edge_space*DX, nx)';
yline = linspace(min(stay) - edge_space*DY, max(stay) + edge_space*DY, ny)';
[xgrid, ygrid] = ndgrid(xline, yline);
[longrid, latgrid] = m_xy2ll(xgrid, ygrid);

fprintf('Inversion grid: %d x %d nodes\n', nx, ny);

for iinv_cell = to_invert
    iinv = iinv_cell{1};

    %% Determine if depth slice or discontinuity parameter
    v_at_depth = isnumeric(iinv);
    if v_at_depth
        param = zcoarse(iinv);
        this_inversion = sprintf('vs%.1f', param);
        rough_scale = rough_scale_base;
        fprintf('\nOn depth slice %d/%d (%.1f km)\n', iinv, nz, param);
    else
        param = iinv;
        this_inversion = sprintf('%s', param);
        rough_scale = rough_scale_params.(char(param));
        fprintf('\nOn parameter %s\n', iinv);
    end

    inv_dir = fullfile(out_dir, this_inversion);
    if ~exist(inv_dir,'dir'), mkdir(inv_dir); end

    fname_surf_vals = fullfile(inv_dir, sprintf('surface_values_V%d.mat', version_surf));
    if exist(fname_surf_vals,'file') && ~re_run
        fprintf('Already have results for %s, skipping\n', this_inversion);
        continue;
    end

    %% Extract station median values for this parameter
    m_simple = zeros(nsta_mdl, 1);
    for imd = 1:nsta_mdl
        mdl = mdls.model{imd};
        if v_at_depth
            [~, iclosest_z] = min(abs(mdl.Z - param));
            m_simple(imd) = mdl.VSav(iclosest_z);
        else
            m_simple(imd) = mdl.([char(iinv) 'av']);
        end
    end

    %% Extract PDFs for this parameter
    if v_at_depth
        iz = iinv;
        pdfs = pdfs_all(1).vs{1};
        for ista = 1:nsta_mdl
            if ~isempty(pdfs_all(ista).vs{iz})
                pdfs(ista) = pdfs_all(ista).vs{iz};
            end
        end
    else
        pdfs = [pdfs_all.(iinv)];
    end

    %% Gaussian interpolation as starting model
    m_interp = griddata(stax, stay, m_simple, xgrid, ygrid, 'cubic');

    mgrid = m_interp;
    mgrid(isnan(mgrid)) = nanmean(mgrid(:));

    for ipt = 1:numel(xgrid)
        distpt = sqrt((stax - xgrid(ipt)).^2 + (stay - ygrid(ipt)).^2);
        distpt = distpt ./ sqrt(DX^2 + DY^2);
        [distpts, distpti] = sort(distpt);
        distpts = distpts(1:nsta_mdl);
        vspt = m_simple(distpti(1:nsta_mdl));
        sta_wt = (1./distpts).^4;
        sta_wt = sta_wt / sum(sta_wt);
        mgrid(ipt) = sum(sta_wt .* vspt);
    end

    a6_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', m_simple, ...
        'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', mgrid, ...
        'fignum', 3, 'title', sprintf('Starting model: %s', this_inversion));
    exportgraphics(gcf, fullfile(inv_dir, sprintf('starting_model_V%d.pdf', version_surf)));

    %% Set up grid roughness
    dx2 = ((xgrid(1:end-2,:) - xgrid(3:end,:))/2).^2;
    dy2 = ((ygrid(:,1:end-2) - ygrid(:,3:end))/2).^2;

    mgrid_start = mgrid;

    %% Prep for efficient inversion
    nmm = 300;
    [pdf_terp, mm_terp, dmm_di] = p_prep_mm_to_pdf(pdfs, nmm);

    [fhand_vec, fhand_mat, grid_terp, nearesti, weighti] = ...
        p_prep_grid_to_sta_interp(xgrid, ygrid, mgrid, stax, stay);

    %% Modify roughness across discontinuities (for velocity inversions only)
    stas_removed = [];
    if v_at_depth
        for this_surf = disconts
            this_surf = this_surf{1};
            disc_file = fullfile(out_dir, this_surf, sprintf('surface_values_V%d.mat', version_surf));
            if ~exist(disc_file, 'file'), continue; end
            zdisc = load(disc_file).mgrid_out;
            gthan = zdisc > param;

            xcomp = (gthan(2:end-1,:) == gthan(3:end,:)) & (gthan(2:end-1,:) == gthan(1:end-2,:));
            ycomp = (gthan(:,2:end-1) == gthan(:,3:end)) & (gthan(:,2:end-1) == gthan(:,1:end-2));
            fprintf('Removing smoothing from %.3f%% x and %.3f%% y cells, %s\n', ...
                100*sum(~xcomp(:))/numel(xcomp), 100*sum(~ycomp(:))/numel(ycomp), this_surf);

            dx2(~xcomp) = inf;
            dy2(~ycomp) = inf;

            for ista = 1:nsta_mdl
                if length(unique(gthan(nearesti(ista,:)))) > 1
                    weighti(ista,:) = 0;
                    grid_terp(:,ista) = 0;
                    stas_removed = [stas_removed; ista];
                end
            end
            fprintf('Removing %d stations due to %s intersection\n', length(stas_removed), this_surf);
        end
    end

    %% Scale roughness by PDF width
    stas_kept = true(nsta_mdl, 1);
    stas_kept(stas_removed) = false;

    bounds_ratio = [0.25, 0.75];
    bounds_stas = zeros(nsta_mdl, 2);
    for ista = 1:nsta_mdl
        mm = pdfs(ista).mm;
        pm = pdfs(ista).pm;
        cpm = cumtrapz(mm, pm);
        get_rid = (cpm > 0.99) | (cpm < 0.01);
        mm = mm(~get_rid); cpm = cpm(~get_rid);
        [~, IA] = unique(cpm);
        cpm = cpm(IA); mm = mm(IA);
        if length(cpm) >= 2
            bounds_stas(ista,:) = interp1(cpm, mm, bounds_ratio, 'spline');
        else
            bounds_stas(ista,:) = [mm(1), mm(end)];
        end
    end
    stas_mod_width = bounds_stas(:,2) - bounds_stas(:,1);
    stas_mod_width(stas_mod_width <= 0) = 1;
    scale_sta_roughness = 1./stas_mod_width;
    mult_roughness = mean(scale_sta_roughness(stas_kept));
    rough_scale = rough_scale * mult_roughness;

    %% Build penalty function handle
    fhand_penalty = @(mgrid)a5_penalty_efficient(mgrid, ...
        pdf_terp, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay, ...
        nearesti, weighti, min(mm_terp), dmm_di, nmm, nsta_mdl);

    pdf_total_max = sum(max(pdf_terp, [], 1));
    fprintf('Max possible PDF sum: %.1f\n', pdf_total_max);

    %% Run fminunc inversion
    opts = optimoptions('fminunc', 'Display', 'iter', ...
        'MaxFunctionEvaluations', inf, 'MaxIterations', max_inv_iterations, ...
        'Algorithm', 'quasi-newton', 'HessianApproximation', 'lbfgs', ...
        'SpecifyObjectiveGradient', false);

    tic;
    ii = 0;
    mgrid_temp = mgrid_start;
    opts_temp = opts;
    ii_vec = 0;
    [pentot_ii, penpdf_ii, penprior_ii, pennorm_ii] = fhand_penalty(mgrid_temp);
    opts_temp.MaxIterations = 0;

    while ii < max_inv_iterations
        if ii < 5
            opts_temp.MaxIterations = 1;
        elseif ii < 15
            opts_temp.MaxIterations = 2;
        elseif ii < 50
            opts_temp.MaxIterations = 6;
        else
            opts_temp.MaxIterations = 15;
        end

        if ii + opts_temp.MaxIterations > max_inv_iterations
            opts_temp.MaxIterations = max(1, max_inv_iterations - ii);
        end

        [mgrid_temp, ~, ~, output_out] = fminunc(fhand_penalty, mgrid_temp, opts_temp);

        ii = ii + output_out.iterations;
        ii_vec(end+1) = ii;
        [pentot_ii(end+1), penpdf_ii(end+1), penprior_ii(end+1), pennorm_ii(end+1)] = ...
            fhand_penalty(mgrid_temp);
    end
    mgrid_out = mgrid_temp;
    toc

    %% Plot inversion progress
    for yscale_type = ["log", "linear"]
        figure(13); clf; hold on; box on;
        set(gca, 'LineWidth', 1.5);
        title(sprintf('Surface inversion: %s (max possible: %.1f)', this_inversion, pdf_total_max), ...
            'FontWeight', 'normal');
        yyaxis('left');
        ylabel('Sum P(m)');
        set(gca, 'YScale', yscale_type, 'YDir', 'reverse');
        grid on;
        plot(ii_vec, -penpdf_ii, '-o');

        yyaxis('right');
        ylabel('Roughness');
        plot(ii_vec, penprior_ii, '-s');
        set(gca, 'YScale', yscale_type);
        xlabel('Iteration');

        exportgraphics(gcf, fullfile(inv_dir, ...
            sprintf('inversion_progress_%s_V%d.pdf', yscale_type, version_surf)));
    end

    %% Plot output surface
    a6_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', m_simple, ...
        'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', mgrid_out, ...
        'fignum', 6, 'title', sprintf('Output: %s', this_inversion));
    exportgraphics(gcf, fullfile(inv_dir, sprintf('surface_output_V%d.pdf', version_surf)));

    %% Save
    save(fullfile(out_dir, 'surface_out_example.mat'), 'longrid', 'latgrid', ...
        'xgrid', 'ygrid', 'llminmax');
    save(fname_surf_vals, 'mgrid_out');
    fprintf('Saved %s\n', fname_surf_vals);
end