clear; close all;

% test_name = 'thick_sed';
% results_dir = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/Results_and_Figures/synthetic_test_correlated_results';
% is_synthetic = true;

test_name = 'D09';
results_dir = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/Results_and_Figures/real_data_correlated_results';
is_synthetic = false;


%% Load chain files
if is_synthetic
    fp = fullfile(results_dir, test_name, sprintf('syn_%s_corr_chain*.mat', test_name));
    files = dir(fp);
    if isempty(files), fp = fullfile(results_dir, sprintf('syn_%s_corr_chain*.mat', test_name)); files = dir(fp); end
else
    fp = fullfile(results_dir, test_name, sprintf('hbi_corr_*%s*.mat', test_name));
    files = dir(fp);
    if isempty(files), fp = fullfile(results_dir, test_name, 'hbi_corr_*.mat'); files = dir(fp); end
end
n_files = length(files);
fprintf('Test: %s (synthetic=%d)\n', test_name, is_synthetic);
fprintf('Found %d chain files\n', n_files);
if n_files == 0, error('No chain files found'); end

all_theta = []; all_sigma = []; all_logL = [];
all_pred_hvsr = []; all_pred_ellip = []; all_pred_cph = []; all_pred_ugr = [];

for k = 1:n_files
    tmp = load(fullfile(files(k).folder, files(k).name));
    r = tmp.results;
    fprintf('  %s: %d samples, best logL=%.2f\n', files(k).name, size(r.samples_theta,1), r.best_logL);
    all_theta = [all_theta; r.samples_theta];
    all_sigma = [all_sigma; r.samples_sigma];
    all_logL  = [all_logL;  r.samples_logL(:)];
    if isfield(r,'samples_pred_hvsr'),  all_pred_hvsr  = [all_pred_hvsr;  r.samples_pred_hvsr];  end
    if isfield(r,'samples_pred_ellip'), all_pred_ellip = [all_pred_ellip; r.samples_pred_ellip]; end
    if isfield(r,'samples_pred_cph'),   all_pred_cph   = [all_pred_cph;   r.samples_pred_cph];   end
    if isfield(r,'samples_pred_ugr'),   all_pred_ugr   = [all_pred_ugr;   r.samples_pred_ugr];   end
    if k == 1
        if isfield(r,'prior'), prior = r.prior; elseif isfield(tmp,'prior'), prior = tmp.prior; end
        if is_synthetic
            if isfield(r,'syn_data_full'), data = r.syn_data_full;
            elseif isfield(r,'syn_data'), data = r.syn_data;
            elseif isfield(tmp,'data'), data = tmp.data; end
            if isfield(r,'true_model'), true_model = r.true_model;
            elseif isfield(tmp,'true_model'), true_model = tmp.true_model; end
        else
            if isfield(tmp,'data'), data = tmp.data; elseif isfield(r,'data'), data = r.data; end
        end
    end
end

theta = all_theta; sigma = all_sigma; logL = all_logL;
Nmods = size(theta, 1);
has_pred_hvsr  = ~isempty(all_pred_hvsr);
has_pred_ellip = ~isempty(all_pred_ellip);
has_pred_cph   = ~isempty(all_pred_cph);
has_pred_ugr   = ~isempty(all_pred_ugr);
pred = struct();
if has_pred_hvsr,  pred.hvsr  = all_pred_hvsr;  end
if has_pred_ellip, pred.ellip = all_pred_ellip; end
if has_pred_cph,   pred.cph   = all_pred_cph;   end
if has_pred_ugr,   pred.ugr   = all_pred_ugr;   end
fprintf('\nCombined: %d chains, %d samples\n', n_files, Nmods);

if is_synthetic && exist('true_model','var') && isfield(true_model,'data_types')
    data_tag = true_model.data_types;
elseif isfield(r,'dtype'), data_tag = r.dtype;
elseif contains(test_name,'hvsr_only'), data_tag = 'hvsr_only';
elseif contains(test_name,'hvsr_ellip'), data_tag = 'hvsr_ellip';
elseif contains(test_name,'hvsr_disp'), data_tag = 'hvsr_disp';
elseif contains(test_name,'ellip_disp'), data_tag = 'ellip_disp';
elseif contains(test_name,'disp_only'), data_tag = 'disp_only';
else, data_tag = 'all';
end


%% Compute Vs profiles
dz = 1/1000;
dvs_brennan = 0.1;
dvs_pdf = 25/1000;
idx_h = prior.idx.h; idx_vs = prior.idx.vs;
h = theta(:,idx_h); vs = theta(:,idx_vs);
zmax = sum(prior.bounds.upper(idx_h));
vsmax = max(prior.bounds.upper(idx_vs));
zmax_actual = max(sum(h,2));
zinterp = 0:dz:zmax; nz = length(zinterp);
L_norm = exp(logL - max(logL));

vs_mod = nan(Nmods,nz); z_mod = nan(Nmods,nz); posterior_mod = zeros(Nmods,nz);
for imod = 1:Nmods
    Igood = theta(imod,:) >= prior.bounds.lower' & theta(imod,:) <= prior.bounds.upper';
    Igood_vs = Igood(idx_vs); Igood_h = [Igood(idx_h), true];
    prior_mod = Igood_h(:) .* Igood_vs(:);
    hmod = [h(imod,:)'; zmax-sum(h(imod,:),2)];
    vsmod = vs(imod,:)';
    modn = [hmod, vsmod*0, vsmod, prior_mod];
    mod_int = layerizemod_interp(modn, zinterp);
    mod_int.prior = mod_int.rho;
    vs_mod(imod,:) = mod_int.vs; z_mod(imod,:) = mod_int.z;
    posterior_mod(imod,:) = L_norm(imod) .* mod_int.prior;
end


%% Brennan-style heatmap (equal weight histogram)
Xvs = 0:dvs_brennan:vsmax*1.05;
hmap = zeros(nz, length(Xvs));
for iz = 1:nz
    hmap(iz,:) = histc(vs_mod(:,iz), Xvs);
end
hmap = hmap / Nmods;
hmap_log = log(hmap);
hmap_log(isinf(hmap_log)) = -20;


%% Likelihood-weighted marginal PDF
vs_edges = 0:dvs_pdf:vsmax*1.05;
vs_vec = 0.5*(vs_edges(1:end-1)+vs_edges(2:end));
marginal_pdf = zeros(nz, length(vs_vec));
for ilay = 1:nz
    ind_bin = discretize(vs_mod(:,ilay), vs_edges);
    marg = zeros(size(vs_vec));
    for ii = 1:Nmods
        if ~isnan(ind_bin(ii))
            marg(ind_bin(ii)) = marg(ind_bin(ii)) + posterior_mod(ii,ilay);
        end
    end
    rs = sum(marg); if rs > 0, marginal_pdf(ilay,:) = marg/rs; end
end

vs_vec_plot = [0, vs_vec];
marginal_pdf_plot = [zeros(nz,1), marginal_pdf];
[xmat, zmat] = meshgrid(vs_vec_plot, zinterp);

valid_rows = sum(marginal_pdf,2) > 0;
med_vs = nan(1,nz); l95_vs = nan(1,nz); u95_vs = nan(1,nz);
l68_vs = nan(1,nz); u68_vs = nan(1,nz);
if any(valid_rows)
    med_vs(valid_rows) = pdf_prctile(marginal_pdf(valid_rows,:), vs_vec, 50);
    l95_vs(valid_rows) = pdf_prctile(marginal_pdf(valid_rows,:), vs_vec, 2.5);
    u95_vs(valid_rows) = pdf_prctile(marginal_pdf(valid_rows,:), vs_vec, 97.5);
    l68_vs(valid_rows) = pdf_prctile(marginal_pdf(valid_rows,:), vs_vec, 16);
    u68_vs(valid_rows) = pdf_prctile(marginal_pdf(valid_rows,:), vs_vec, 84);
end

if is_synthetic && exist('true_model','var')
    h_true = true_model.theta(idx_h); vs_true = true_model.theta(idx_vs);
    z_true = [0; cumsum(h_true)];
end
[~, ibest] = max(logL);
vs_map = vs_mod(ibest,:);

T_hvsr=[]; T_ellip=[]; T_ugr=[]; T_cph=[];
if isfield(data,'hvsr') && isfield(data.hvsr,'f'), T_hvsr = 1./data.hvsr.f(:); end
if isfield(data,'ellip') && isfield(data.ellip,'T'), T_ellip = data.ellip.T(:); end
if isfield(data,'ugr') && isfield(data.ugr,'T'), T_ugr = data.ugr.T(:); end
if isfield(data,'cph') && isfield(data.cph,'T'), T_cph = data.cph.T(:); end


%% Color setup
clrbase = viridis;
cmap = [flipud([linspace(clrbase(1,1),0.8,10)',linspace(clrbase(1,2),0.8,10)',linspace(clrbase(1,3),0.8,10)']); clrbase];
clims = [prctile(logL,25) max(logL)];
if clims(1)>=clims(2) || any(~isfinite(clims)), clims = [max(logL)-10, max(logL)]; end
cvals = linspace(clims(1),clims(2),size(cmap,1));
[~,isrt] = sort(logL);

clims_post = [-4 0];
log10_pdf = log10(marginal_pdf_plot);
log10_pdf(~isfinite(log10_pdf)) = clims_post(1);

plot_true = is_synthetic && exist('vs_true','var');
logL_thresh = prctile(logL,25);

savedir = fullfile(results_dir, test_name, 'test final model overlay figures');
if ~exist(savedir,'dir'), mkdir(savedir); end

zinterp_col = zinterp(:);
med_col = med_vs(:); l95_col = l95_vs(:); u95_col = u95_vs(:);
l68_col = l68_vs(:); u68_col = u68_vs(:);


%% Interface PPD
disc_dz = 0.05;
disc_z_edges = 0:disc_dz:zmax;
disc_z_vec = 0.5*(disc_z_edges(1:end-1)+disc_z_edges(2:end));
disc_pdf = zeros(size(disc_z_vec));
for imod = 1:Nmods
    if logL(imod) < logL_thresh, continue; end
    ifd = cumsum(h(imod,:)); vsi = vs(imod,:);
    for id = 1:length(ifd)
        dv = abs(vsi(id+1)-vsi(id));
        if dv < 0.1, continue; end
        ib = discretize(ifd(id), disc_z_edges);
        if ~isempty(ib) && ~isnan(ib)
            disc_pdf(ib) = disc_pdf(ib) + L_norm(imod)*dv;
        end
    end
end
ds = sum(disc_pdf); if ds>0, disc_pdf = disc_pdf/ds; end
if plot_true, z_true_ifaces = cumsum(h_true); end

zoom_depths = [0.5, 1, 5]; zoom_vs_max = [3, 3.5, 4];
zoom_titles = {'0-500m','0-1km','0-5km'};


%% ====== SPAGHETTI FIGURES (1001, 1003-1005) ======

%% Figure 1001: Spaghetti Vs + PPD (full depth)
figure(1001); clf;
set(gcf,'Position',[47 430 900 700],'Color','w','Visible','off'); colormap(cmap);
axV = axes('Position',[0.10 0.08 0.48 0.85]); box on; hold on;
for imod = isrt(:)'
    if logL(imod)<logL_thresh, continue; end
    [~,ci]=min(abs(cvals-logL(imod)));
    plot(vs_mod(imod,:),z_mod(imod,:),'-','color',cmap(ci,:),'linewidth',3);
end
if plot_true, plot_true_model(vs_true,z_true); end
set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse','Layer','top');
pp=get(gca,'Position'); cbb=colorbar; set(gca,'Position',pp);
ylabel(cbb,'log(Likelihood)'); set(cbb,'linewidth',1.5,'fontsize',14);
xlim([0 max(vs_vec)]); ylim([0 zmax_actual]); caxis(clims);
xlabel('Vs (km/s)'); ylabel('Depth (km)');
title(sprintf('Vs Profiles [%s]',data_tag),'FontSize',14);

axD = axes('Position',[0.68 0.08 0.25 0.85]); box on; hold on;
barh(disc_z_vec,disc_pdf,1.0,'FaceColor',[0.9 0 0],'EdgeColor','none');
if plot_true, plot_true_ifaces(z_true_ifaces, zmax_actual); end
set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse','Layer','top');
ylim([0 zmax_actual]); xlabel('Interface PPD'); set(gca,'YTickLabel',[]);
linkaxes([axV,axD],'y');

%% Figures 1003-1005: Zoomed spaghetti + PPD
for iz = 1:3
    figure(1002+iz); clf;
    set(gcf,'Position',[50 400 900 700],'Color','w','Visible','off'); colormap(cmap);
    axA=axes('Position',[0.10 0.08 0.48 0.85]); box on; hold on;
    for imod=isrt(:)'
        if logL(imod)<logL_thresh, continue; end
        [~,ci]=min(abs(cvals-logL(imod)));
        plot(vs_mod(imod,:),z_mod(imod,:),'-','color',cmap(ci,:),'linewidth',3);
    end
    if plot_true, plot_true_model(vs_true,z_true); end
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse','Layer','top');
    pp=get(gca,'Position'); cbb=colorbar; set(gca,'Position',pp);
    ylabel(cbb,'log(Likelihood)'); set(cbb,'linewidth',1.5,'fontsize',14);
    xlim([0 zoom_vs_max(iz)]); ylim([0 zoom_depths(iz)]); caxis(clims);
    xlabel('Vs (km/s)'); ylabel('Depth (km)');
    title(sprintf('Vs Profiles (%s) [%s]',zoom_titles{iz},data_tag),'FontSize',14);

    axB=axes('Position',[0.68 0.08 0.25 0.85]); box on; hold on;
    barh(disc_z_vec,disc_pdf,1.0,'FaceColor',[0.9 0 0],'EdgeColor','none');
    if plot_true, plot_true_ifaces(z_true_ifaces, zoom_depths(iz)); end
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse','Layer','top');
    ylim([0 zoom_depths(iz)]); xlabel('Interface PPD'); set(gca,'YTickLabel',[]);
    linkaxes([axA,axB],'y');
end


%% ====== BRENNAN HEATMAP FIGURES (1015-1020) ======

%% Build NaN-safe fill vectors for credible intervals on heatmap
iv_ok = ~isnan(l68_col) & ~isnan(u68_col);
fill_x68 = [l68_col(iv_ok); flipud(u68_col(iv_ok))];
fill_z68 = [zinterp_col(iv_ok); flipud(zinterp_col(iv_ok))];
iv_ok95 = ~isnan(l95_col) & ~isnan(u95_col);
fill_x95 = [l95_col(iv_ok95); flipud(u95_col(iv_ok95))];
fill_z95 = [zinterp_col(iv_ok95); flipud(zinterp_col(iv_ok95))];

%% Figure 1015: Brennan heatmap + final model + PPD (full depth)
figure(1015); clf;
set(gcf,'Position',[50 100 900 800],'Color','w','Visible','off');
axV = axes('Position',[0.10 0.08 0.48 0.85],'linewidth',1.5); box on; hold on;
contourf(axV, Xvs, zinterp, hmap_log, [-5:0.1:-0.1], 'EdgeColor','none');
fill(axV, fill_x68, fill_z68, [0.7 0.7 0.7], 'LineWidth', 1.5,...
    'EdgeColor', [0.6 0.6 0.6], 'FaceAlpha', 0.4);
plot(axV, l95_vs, zinterp, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
plot(axV, u95_vs, zinterp, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
plot(axV, med_vs, zinterp, '--k', 'LineWidth', 2);
if plot_true, plot_true_model(vs_true,z_true); end
set(axV,'ydir','reverse','FontSize',15,'linewidth',1.5,'Color','none');
xlim(axV,[0 vsmax*1.05]); ylim(axV,[0 zmax_actual]);
xlabel(axV,'Vs (km/s)','FontSize',18); ylabel(axV,'Depth (km)','FontSize',18);
caxis(axV,[-5 0]);
hcb = colorbar(axV);
set(axV,'Position',[0.10 0.08 0.48 0.85]);
set(hcb,'Position',[0.59 0.08 0.02 0.85],'linewidth',1.5,'fontsize',14);
hcby = ylabel(hcb,'log(Probability)','FontWeight','bold','FontSize',14);
set(hcby,'rotation',270,'VerticalAlignment','bottom');
title(axV, sprintf('Posterior Vs [%s]',data_tag),'FontSize',18,'FontWeight','bold');

axD = axes('Position',[0.68 0.08 0.25 0.85],'linewidth',1.5); box on; hold on;
barh(disc_z_vec,disc_pdf,1.0,'FaceColor',[0.9 0 0],'EdgeColor','none');
if plot_true, plot_true_ifaces(z_true_ifaces, zmax_actual); end
set(axD,'ydir','reverse','FontSize',15,'linewidth',1.5);
ylim(axD,[0 zmax_actual]); xlabel(axD,'Interface PPD'); set(axD,'YTickLabel',[]);
linkaxes([axV,axD],'y');

%% Figures 1016-1018: Zoomed Brennan heatmap + final model + PPD
for iz = 1:3
    figure(1015+iz); clf;
    set(gcf,'Position',[50 100 900 800],'Color','w','Visible','off');
    axV = axes('Position',[0.10 0.08 0.48 0.85],'linewidth',1.5); box on; hold on;
    contourf(axV, Xvs, zinterp, hmap_log, [-5:0.1:-0.1], 'EdgeColor','none');
    fill(axV, fill_x68, fill_z68, [0.7 0.7 0.7], 'LineWidth', 1.5,...
        'EdgeColor', [0.6 0.6 0.6], 'FaceAlpha', 0.4);
    plot(axV, l95_vs, zinterp, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    plot(axV, u95_vs, zinterp, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    plot(axV, med_vs, zinterp, '--k', 'LineWidth', 2);
    if plot_true, plot_true_model(vs_true,z_true); end
    set(axV,'ydir','reverse','FontSize',15,'linewidth',1.5,'Color','none');
    xlim(axV,[0 zoom_vs_max(iz)]); ylim(axV,[0 zoom_depths(iz)]);
    xlabel(axV,'Vs (km/s)','FontSize',18); ylabel(axV,'Depth (km)','FontSize',18);
    caxis(axV,[-5 0]);
    hcb = colorbar(axV);
    set(axV,'Position',[0.10 0.08 0.48 0.85]);
    set(hcb,'Position',[0.59 0.08 0.02 0.85],'linewidth',1.5,'fontsize',14);
    hcby = ylabel(hcb,'log(Probability)','FontWeight','bold','FontSize',14);
    set(hcby,'rotation',270,'VerticalAlignment','bottom');
    title(axV, sprintf('Posterior Vs (%s) [%s]',zoom_titles{iz},data_tag),'FontSize',18,'FontWeight','bold');

    axD = axes('Position',[0.68 0.08 0.25 0.85],'linewidth',1.5); box on; hold on;
    barh(disc_z_vec,disc_pdf,1.0,'FaceColor',[0.9 0 0],'EdgeColor','none');
    if plot_true, plot_true_ifaces(z_true_ifaces, zoom_depths(iz)); end
    set(axD,'ydir','reverse','FontSize',15,'linewidth',1.5);
    ylim(axD,[0 zoom_depths(iz)]); xlabel(axD,'Interface PPD'); set(axD,'YTickLabel',[]);
    linkaxes([axV,axD],'y');
end


%% ====== POSTERIOR PDF FIGURES (1002, 1006-1008) ======

%% Figure 1002: Full depth posterior PDF + CI
figure(1002); clf;
set(gcf,'Position',[1 425 950 527],'Color','w','Visible','off'); colormap(cmap);
axes('Position',[0.07 0.12 0.33 0.76]); box on; hold on;
surface(xmat,zmat,zeros(size(zmat)),log10_pdf,'LineStyle','none'); shading interp;
if plot_true, plot_true_model(vs_true,z_true); end
plot(med_vs,zinterp,'-w','LineWidth',2);
plot(l68_vs,zinterp,'--w','LineWidth',1.5); plot(u68_vs,zinterp,'--w','LineWidth',1.5);
set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse','Layer','top');
cb=colorbar; set(gca,'Position',[0.07 0.12 0.33 0.76]); set(cb,'Position',[0.41 0.12 0.02 0.76]);
ylabel(cb,'log_{10}(Probability)'); set(cb,'linewidth',1.5,'fontsize',14);
xlim([0 max(vs_vec)]); ylim([0 zmax_actual]); caxis(clims_post);
xlabel('Vs (km/s)'); ylabel('Depth (km)'); title('Posterior Probability','FontSize',14);

axes('Position',[0.56 0.12 0.38 0.76]); box on; hold on;
h95=fill([l95_col;flipud(u95_col)],[zinterp_col;flipud(zinterp_col)],[0.8 0.8 1.0],'EdgeColor','none','FaceAlpha',0.4);
h68=fill([l68_col;flipud(u68_col)],[zinterp_col;flipud(zinterp_col)],[0.6 0.6 1.0],'EdgeColor','none','FaceAlpha',0.5);
hmed=plot(med_col,zinterp_col,'-b','LineWidth',2);
if plot_true, htrue=plot_true_model(vs_true,z_true); end
set(gca,'FontSize',16,'linewidth',1.5,'ydir','reverse','Layer','top');
xlim([0 max(vs_vec)]); ylim([0 zmax_actual]);
xlabel('Vs (km/s)'); ylabel('Depth (km)'); title('Median + Confidence','FontSize',16);
if plot_true, legend([h95,h68,hmed,htrue],'95% CI','68% CI','Median','True','Location','southwest');
else, legend([h95,h68,hmed],'95% CI','68% CI','Median','Location','southwest'); end
sgtitle(sprintf('Posterior [%s]: %s',data_tag,test_name),'FontSize',16,'FontWeight','bold');

%% Figures 1006-1008: Zoomed posterior PDF + CI
for iz = 1:3
    figure(1005+iz); clf;
    set(gcf,'Position',[100 400 950 500],'Color','w','Visible','off'); colormap(cmap);
    axes('Position',[0.07 0.12 0.33 0.76]); box on; hold on;
    surface(xmat,zmat,zeros(size(zmat)),log10_pdf,'LineStyle','none'); shading interp;
    if plot_true, plot_true_model(vs_true,z_true); end
    plot(med_vs,zinterp,'-w','LineWidth',2);
    plot(l68_vs,zinterp,'--w','LineWidth',1.5); plot(u68_vs,zinterp,'--w','LineWidth',1.5);
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse','Layer','top');
    cb=colorbar; set(gca,'Position',[0.07 0.12 0.33 0.76]); set(cb,'Position',[0.41 0.12 0.02 0.76]);
    ylabel(cb,'log_{10}(Probability)'); set(cb,'linewidth',1.5,'fontsize',14);
    xlim([0 zoom_vs_max(iz)]); ylim([0 zoom_depths(iz)]); caxis(clims_post);
    xlabel('Vs (km/s)'); ylabel('Depth (km)'); title(sprintf('Posterior (%s)',zoom_titles{iz}),'FontSize',14);

    axes('Position',[0.56 0.12 0.38 0.76]); box on; hold on;
    h95z=fill([l95_col;flipud(u95_col)],[zinterp_col;flipud(zinterp_col)],[0.8 0.8 1.0],'EdgeColor','none','FaceAlpha',0.4);
    h68z=fill([l68_col;flipud(u68_col)],[zinterp_col;flipud(zinterp_col)],[0.6 0.6 1.0],'EdgeColor','none','FaceAlpha',0.5);
    hmedz=plot(med_col,zinterp_col,'-b','LineWidth',2);
    if plot_true, htruez=plot_true_model(vs_true,z_true); end
    set(gca,'FontSize',16,'linewidth',1.5,'ydir','reverse','Layer','top');
    xlim([0 zoom_vs_max(iz)]); ylim([0 zoom_depths(iz)]);
    xlabel('Vs (km/s)'); ylabel('Depth (km)'); title(sprintf('Median + CI (%s)',zoom_titles{iz}),'FontSize',16);
    if plot_true, legend([h95z,h68z,hmedz,htruez],'95% CI','68% CI','Median','True','Location','southeast');
    else, legend([h95z,h68z,hmedz],'95% CI','68% CI','Median','Location','southeast'); end
end


%% ====== OBSERVED vs PREDICTED (1010) ======
figure(1010); clf;
set(gcf,'Position',[50 300 900 700],'Color','w','Visible','off'); colormap(cmap);
ipanel = 0;
n_errbar = 50;

if has_pred_hvsr && ~isempty(T_hvsr)
    ipanel=ipanel+1; subplot(2,2,ipanel); box on; hold on;
    for imod=isrt(:)', if logL(imod)<logL_thresh, continue; end
        [~,ci]=min(abs(cvals-logL(imod)));
        plot(T_hvsr,pred.hvsr(imod,:)','-','color',[cmap(ci,:) 0.3],'linewidth',1);
    end
    plot(T_hvsr,pred.hvsr(ibest,:)','-k','LineWidth',2);
    if plot_true && isfield(true_model,'clean') && isfield(true_model.clean,'hvsr')
        plot(T_hvsr,true_model.clean.hvsr,'r-','LineWidth',2); end
    eidx = round(linspace(1,length(T_hvsr),min(n_errbar,length(T_hvsr))));
    errorbar(T_hvsr(eidx),data.hvsr.obs(eidx),data.hvsr.sig(eidx),'.k','linewidth',2);
    xlabel('Period (s)'); ylabel('HVSR');
    set(gca,'FontSize',14,'linewidth',1.5,'Layer','top','XScale','log'); title('HVSR','FontSize',14);
end

if has_pred_ellip && ~isempty(T_ellip)
    ipanel=ipanel+1; subplot(2,2,ipanel); box on; hold on;
    for imod=isrt(:)', if logL(imod)<logL_thresh, continue; end
        [~,ci]=min(abs(cvals-logL(imod)));
        plot(T_ellip,pred.ellip(imod,:)','-','color',[cmap(ci,:) 0.3],'linewidth',1);
    end
    plot(T_ellip,pred.ellip(ibest,:)','-k','LineWidth',2);
    if plot_true && isfield(true_model,'clean') && isfield(true_model.clean,'ellip')
        plot(T_ellip,true_model.clean.ellip,'r-','LineWidth',2); end
    eidx = round(linspace(1,length(T_ellip),min(n_errbar,length(T_ellip))));
    errorbar(T_ellip(eidx),data.ellip.obs(eidx),data.ellip.sig(eidx),'.k','linewidth',2);
    xlabel('Period (s)'); ylabel('Ellipticity');
    set(gca,'FontSize',14,'linewidth',1.5,'Layer','top'); title('Rayleigh Wave Ellipticity','FontSize',14);
end

if has_pred_cph && ~isempty(T_cph)
    ipanel=ipanel+1; subplot(2,2,ipanel); box on; hold on;
    for imod=isrt(:)', if logL(imod)<logL_thresh, continue; end
        [~,ci]=min(abs(cvals-logL(imod)));
        plot(T_cph,pred.cph(imod,:)','-','color',[cmap(ci,:) 0.3],'linewidth',1);
    end
    plot(T_cph,pred.cph(ibest,:)','-k','LineWidth',2);
    if plot_true && isfield(true_model,'clean') && isfield(true_model.clean,'cph')
        plot(T_cph,true_model.clean.cph,'r-','LineWidth',2); end
    errorbar(T_cph,data.cph.obs,data.cph.sig,'.k','linewidth',2);
    xlabel('Period (s)'); ylabel('Phase Velocity (km/s)');
    set(gca,'FontSize',14,'linewidth',1.5,'Layer','top'); title('Phase Velocity Dispersion','FontSize',14);
end

if has_pred_ugr && ~isempty(T_ugr)
    ipanel=ipanel+1; subplot(2,2,ipanel); box on; hold on;
    for imod=isrt(:)', if logL(imod)<logL_thresh, continue; end
        [~,ci]=min(abs(cvals-logL(imod)));
        plot(T_ugr,pred.ugr(imod,:)','-','color',[cmap(ci,:) 0.3],'linewidth',1);
    end
    plot(T_ugr,pred.ugr(ibest,:)','-k','LineWidth',2);
    if plot_true && isfield(true_model,'clean') && isfield(true_model.clean,'ugr')
        plot(T_ugr,true_model.clean.ugr,'r-','LineWidth',2); end
    errorbar(T_ugr,data.ugr.obs,data.ugr.sig,'.k','linewidth',2);
    xlabel('Period (s)'); ylabel('Group Velocity (km/s)');
    set(gca,'FontSize',14,'linewidth',1.5,'Layer','top'); title('Group Velocity Dispersion','FontSize',14);
end
sgtitle(sprintf('Observed vs Predicted [%s]: %s',data_tag,test_name),'FontSize',15,'FontWeight','bold');


%% ====== SAVE ALL ======
fig_list = {
    1001, 'Vs_spaghetti_full_depth';
    1002, 'Posterior_PDF_full_depth';
    1003, 'Vs_spaghetti_0-500m';
    1004, 'Vs_spaghetti_0-1km';
    1005, 'Vs_spaghetti_0-5km';
    1006, 'Posterior_PDF_0-500m';
    1007, 'Posterior_PDF_0-1km';
    1008, 'Posterior_PDF_0-5km';
    1010, 'Observed_vs_Predicted';
    1015, 'Brennan_heatmap_full_depth';
    1016, 'Brennan_heatmap_0-500m';
    1017, 'Brennan_heatmap_0-1km';
    1018, 'Brennan_heatmap_0-5km';
};

for ii = 1:size(fig_list,1)
    fn = fig_list{ii,1}; nm = sprintf('%s_%s', fig_list{ii,2}, data_tag);
    f = figure(fn); drawnow;
    exportgraphics(f, fullfile(savedir,[nm '.png']), 'Resolution', 300);
    exportgraphics(f, fullfile(savedir,[nm '.pdf']), 'ContentType','image','Resolution',300);
    fprintf('Saved %d: %s\n', fn, nm); close(f);
end
fprintf('All saved to: %s\n', savedir);


%% ====== LOCAL FUNCTIONS ======
function h_out = plot_true_model(vs_true, z_true)
    h_out = [];
    for i = 1:length(vs_true)
        zt = z_true(i);
        if i < length(vs_true), zb = z_true(i+1); else, zb = z_true(end)+5; end
        hh = plot([vs_true(i) vs_true(i)],[zt zb],'r-','LineWidth',2.5);
        if isempty(h_out), h_out = hh; end
        if i < length(vs_true), plot([vs_true(i) vs_true(i+1)],[zb zb],'r-','LineWidth',2.5); end
    end
end

function plot_true_ifaces(z_ifaces, max_depth)
    for it = 1:length(z_ifaces)
        if z_ifaces(it) <= max_depth
            plot(xlim, [z_ifaces(it) z_ifaces(it)], 'b-', 'LineWidth', 1.5);
        end
    end
end