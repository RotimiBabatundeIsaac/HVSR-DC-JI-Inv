function [log_L, log_L_components, pred, Phi, N_d] = c3_compute_likelihood(theta, sigma_e, data, prior)

persistent paths_set cps_bin_set hvf_bin_set

cps_bin = prior.paths.cps_bin;
hvf_bin = prior.paths.hv_dfa_bin;
func_path = fullfile(prior.paths.base, 'AGU_models', 'functions');

need_reset = isempty(paths_set) || ...
             ~strcmp(cps_bin_set, cps_bin) || ...
             ~strcmp(hvf_bin_set, hvf_bin);

if need_reset
    new_path = [cps_bin, ':', hvf_bin, ':/usr/bin:/bin:/usr/sbin:/sbin'];
    setenv('PATH', new_path);
    if ~contains(path, func_path)
        addpath(func_path);
    end
    cps_bin_set = cps_bin;
    hvf_bin_set = hvf_bin;
    paths_set = true;
end

h  = theta(prior.idx.h);
vs = theta(prior.idx.vs);
model_matrix = build_model_matrix(h, vs);

config   = build_forward_config(data);
fwd_data = build_forward_data(data);
pred_raw = forward_model(model_matrix, fwd_data, config);

log_L_components = zeros(4, 1);
Phi  = NaN(4, 1);
N_d  = zeros(4, 1);
pred = struct('hvsr', [], 'ellip', [], 'cph', [], 'ugr', []);

if config.use_hvsr
    [log_L_components(1), Phi(1), N_d(1), pred.hvsr] = compute_likelihood_partial_prefix(...
        data, pred_raw, 'hvsr', 'f', sigma_e(1), 0.25);
else
    log_L_components(1) = 0;
end

if config.use_ellip
    [log_L_components(2), Phi(2), N_d(2), pred.ellip] = compute_likelihood_partial_prefix(...
        data, pred_raw, 'ellip', 'T', sigma_e(2), 0.50);
else
    log_L_components(2) = 0;
end

if config.use_cph
    [log_L_components(3), Phi(3), N_d(3), pred.cph] = compute_likelihood_component(...
        data, pred_raw, 'cph', 'T', sigma_e(3), 0.50);
else
    log_L_components(3) = 0;
end

if config.use_ugr
    [log_L_components(4), Phi(4), N_d(4), pred.ugr] = compute_likelihood_component(...
        data, pred_raw, 'ugr', 'T', sigma_e(4), 0.50);
else
    log_L_components(4) = 0;
end

log_L = sum(log_L_components);
if ~isfinite(log_L)
    log_L = -1e10;
end

end


function [log_L, Phi, N_valid, pred_out] = compute_likelihood_partial_prefix(...
    data, pred_raw, field, x_field, sigma_d, min_fraction)

pred_out = [];
log_L    = -1e10;
Phi      = Inf;
N_valid  = 0;

if ~isfield(data, field), return; end
if ~isfield(data.(field), x_field) || isempty(data.(field).(x_field)), return; end
if ~isfield(pred_raw, field), return; end
if ~isfield(pred_raw.(field), 'status') || pred_raw.(field).status ~= 0, return; end
if ~isfield(pred_raw.(field), 'val'), return; end

pred_out = pred_raw.(field).val(:);
obs      = data.(field).obs(:);
sig      = data.(field).sig(:);

n_pred = numel(pred_out);
n_obs  = numel(obs);
n_sig  = numel(sig);

if n_obs == 0 || n_pred == 0 || n_sig == 0
    pred_out = [];
    return;
end

if n_pred ~= n_obs || n_sig ~= n_obs
    pred_out = [];
    return;
end

valid_obs  = isfinite(obs) & isfinite(sig) & (sig > 0);
N_obs_good = sum(valid_obs);
if N_obs_good == 0
    pred_out = [];
    return;
end

finite_at_valid = valid_obs & isfinite(pred_out);
if ~any(finite_at_valid)
    pred_out = [];
    return;
end
last_finite_idx = find(finite_at_valid, 1, 'last');

keep = valid_obs;
if last_finite_idx < n_obs
    keep((last_finite_idx+1):end) = false;
end

if any(~isfinite(pred_out(keep)))
    pred_out = [];
    return;
end

N_valid = sum(keep);
if N_valid == 0
    pred_out = [];
    return;
end

if (N_valid / N_obs_good) < min_fraction
    pred_out = [];
    log_L = -1e10;
    Phi = Inf;
    N_valid = 0;
    return;
end

obs_v  = obs(keep);
pred_v = pred_out(keep);
sig_v  = sig(keep);

r_raw = (obs_v - pred_v) ./ sig_v;
Phi   = sum(r_raw.^2);

sigma_d   = max(sigma_d, 1e-12);
sigma_eff = sigma_d .* sig_v;

log_L = -0.5 * N_valid * log(2*pi) ...
        - sum(log(sigma_eff)) ...
        - 0.5 * (Phi / (sigma_d^2));

if ~isfinite(log_L)
    log_L = -1e10;
end

end


function [log_L, Phi, N_valid, pred_out] = compute_likelihood_component(...
    data, pred_raw, field, x_field, sigma_d, min_fraction)

pred_out = [];
log_L = -1e10;
Phi = Inf;
N_valid = 0;

if ~isfield(data, field), return; end
if ~isfield(data.(field), x_field) || isempty(data.(field).(x_field)), return; end
if ~isfield(pred_raw, field), return; end
if ~isfield(pred_raw.(field), 'status') || pred_raw.(field).status ~= 0, return; end
if ~isfield(pred_raw.(field), 'val') || isempty(pred_raw.(field).val), return; end

obs = data.(field).obs(:);
sig = data.(field).sig(:);
n_obs = numel(obs);

pred_val = pred_raw.(field).val(:);
n_pred = numel(pred_val);

pred_out = NaN(n_obs, 1);
n_use = min(n_pred, n_obs);
pred_out(1:n_use) = pred_val(1:n_use);

valid_obs = isfinite(obs) & isfinite(sig) & (sig > 0);
N_obs_good = sum(valid_obs);
if N_obs_good == 0
    pred_out = [];
    return;
end

computable = valid_obs & isfinite(pred_out);
N_valid = sum(computable);

if N_valid == 0
    pred_out = [];
    return;
end

if (N_valid / N_obs_good) < min_fraction
    pred_out = [];
    log_L = -1e10;
    Phi = Inf;
    N_valid = 0;
    return;
end

obs_v = obs(computable);
pred_v = pred_out(computable);
sig_v = sig(computable);

r_raw = (obs_v - pred_v) ./ sig_v;
Phi = sum(r_raw.^2);

sigma_d = max(sigma_d, 1e-12);
sigma_eff = sigma_d .* sig_v;

log_L = -0.5 * N_valid * log(2*pi) ...
        - sum(log(sigma_eff)) ...
        - 0.5 * (Phi / (sigma_d^2));

if ~isfinite(log_L)
    log_L = -1e10;
end

end


function model_matrix = build_model_matrix(h, vs)

n_layers = numel(vs);
n_thick = numel(h);
model_matrix = zeros(n_layers, 4);

for i = 1:n_layers
    if i <= n_thick
        model_matrix(i, 1) = h(i);
    else
        model_matrix(i, 1) = 0;
    end
    model_matrix(i, 3) = vs(i);
    [model_matrix(i, 2), model_matrix(i, 4)] = brocher2005(vs(i));
end

end


function [vp, rho] = brocher2005(vs)

vp  = 0.9409 + 2.0947*vs - 0.8206*vs^2 + 0.2683*vs^3 - 0.0251*vs^4;
rho = 1.6612*vp - 0.4721*vp^2 + 0.0671*vp^3 - 0.0043*vp^4 + 0.000106*vp^5;
vp  = max(vp, vs * 1.5);
rho = max(rho, 1.5);

end


function config = build_forward_config(data)

config = struct();
config.use_hvsr  = isfield(data, 'hvsr')  && isfield(data.hvsr,  'f') && ~isempty(data.hvsr.f);
config.use_ellip = isfield(data, 'ellip') && isfield(data.ellip, 'T') && ~isempty(data.ellip.T);
config.use_cph   = isfield(data, 'cph')   && isfield(data.cph,   'T') && ~isempty(data.cph.T);
config.use_ugr   = isfield(data, 'ugr')   && isfield(data.ugr,   'T') && ~isempty(data.ugr.T);

config.hvsr.nf   = 100;
config.hvsr.fmin = 1/3;
config.hvsr.fmax = 25;
config.hvsr.nmr  = 10;
config.hvsr.nml  = 10;
config.hvsr.nks  = 1000;

config.ellip.nf   = 100;
config.disp.nmode = 0;

end


function fwd_data = build_forward_data(data)

fwd_data = struct();

if isfield(data, 'hvsr') && isfield(data.hvsr, 'f') && ~isempty(data.hvsr.f)
    fwd_data.hvsr.f   = data.hvsr.f(:);
    fwd_data.hvsr.obs = data.hvsr.obs(:);
    fwd_data.hvsr.sig = data.hvsr.sig(:);
end

if isfield(data, 'ellip') && isfield(data.ellip, 'T') && ~isempty(data.ellip.T)
    fwd_data.ellip.T   = data.ellip.T(:);
    fwd_data.ellip.obs = data.ellip.obs(:);
    fwd_data.ellip.sig = data.ellip.sig(:);
end

if isfield(data, 'cph') && isfield(data.cph, 'T') && ~isempty(data.cph.T)
    fwd_data.cph.T   = data.cph.T(:);
    fwd_data.cph.obs = data.cph.obs(:);
    fwd_data.cph.sig = data.cph.sig(:);
end

if isfield(data, 'ugr') && isfield(data.ugr, 'T') && ~isempty(data.ugr.T)
    fwd_data.ugr.T   = data.ugr.T(:);
    fwd_data.ugr.obs = data.ugr.obs(:);
    fwd_data.ugr.sig = data.ugr.sig(:);
end

end