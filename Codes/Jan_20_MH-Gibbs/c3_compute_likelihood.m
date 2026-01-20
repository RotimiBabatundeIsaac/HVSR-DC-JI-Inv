function [log_L, log_L_components, pred] = c3_compute_likelihood(theta, sigma_e, data, prior)
% Compute log-likelihood for all four data types
%
% IMPORTANT (consistency with updated c8):
%   - data.*.sig are per-point observational uncertainties
%   - sigma_e(d) is a unitless noise multiplier for each datatype
%   - effective std for point i is: sigma_eff(i) = sigma_e(d) * data.*.sig(i)

persistent paths_set
if isempty(paths_set)
    cps_bin = prior.paths.cps_bin;
    hvf_bin = prior.paths.hv_dfa_bin;
    current_path = getenv('PATH');
    if ~contains(current_path, cps_bin)
        setenv('PATH', [hvf_bin, ':', cps_bin, ':', current_path]);
    end
    paths_set = true;
end

h  = theta(1:5);
vs = theta(6:11);
model = build_model_matrix(h, vs);

config   = build_forward_config(data);
fwd_data = build_forward_data(data);
pred_raw = forward_model(model, fwd_data, config);

log_L_components = zeros(4, 1);
pred = struct();

% -------------------- 1) HVSR --------------------
if isfield(data, 'hvsr') && ~isempty(data.hvsr.f) && pred_raw.hvsr.status == 0
    pred.hvsr = pred_raw.hvsr.val(:);

    % Diagnostic: where prediction becomes non-finite
    try
        f_obs = data.hvsr.f(:);
        ok = isfinite(f_obs) & isfinite(pred.hvsr);
        if any(ok)
            fprintf('[c3] HVSR finite: %d/%d, max finite f = %.2f Hz\n', sum(ok), numel(pred.hvsr), max(f_obs(ok)));
        else
            fprintf('[c3] HVSR finite: 0/%d (all NaN/Inf)\n', numel(pred.hvsr));
        end
    catch
    end

    obs = data.hvsr.obs(:);
    sig = data.hvsr.sig(:);

    valid = isfinite(obs) & isfinite(sig) & (sig > 0) & isfinite(pred.hvsr);
    n_valid = sum(valid);
    n_total = numel(obs);

    if n_valid > 0
        log_L_components(1) = gaussian_log_likelihood_hetero(obs(valid), pred.hvsr(valid), sig(valid), sigma_e(1));
        if n_valid < n_total
            log_L_components(1) = log_L_components(1) * (n_valid / n_total);
        end
    else
        log_L_components(1) = -1e10;
    end
else
    log_L_components(1) = -1e10;
    pred.hvsr = [];
end

% -------------------- 2) Ellipticity --------------------
if isfield(data, 'ellip') && ~isempty(data.ellip.T) && pred_raw.ellip.status == 0
    pred.ellip = pred_raw.ellip.val(:);

    obs = data.ellip.obs(:);
    sig = data.ellip.sig(:);

    valid = isfinite(obs) & isfinite(sig) & (sig > 0) & isfinite(pred.ellip);
    n_valid = sum(valid);
    n_total = numel(obs);

    if n_valid > 0
        log_L_components(2) = gaussian_log_likelihood_hetero(obs(valid), pred.ellip(valid), sig(valid), sigma_e(2));
        if n_valid < n_total
            log_L_components(2) = log_L_components(2) * (n_valid / n_total);
        end
    else
        log_L_components(2) = -1e10;
    end
else
    log_L_components(2) = -1e10;
    pred.ellip = [];
end

% -------------------- 3) Phase velocity --------------------
if isfield(data, 'cph') && ~isempty(data.cph.T) && pred_raw.cph.status == 0
    pred.cph = pred_raw.cph.val(:);

    obs = data.cph.obs(:);
    sig = data.cph.sig(:);

    valid = isfinite(obs) & isfinite(sig) & (sig > 0) & isfinite(pred.cph);
    n_valid = sum(valid);
    n_total = numel(obs);

    if n_valid > 0
        log_L_components(3) = gaussian_log_likelihood_hetero(obs(valid), pred.cph(valid), sig(valid), sigma_e(3));
        if n_valid < n_total
            log_L_components(3) = log_L_components(3) * (n_valid / n_total);
        end
    else
        log_L_components(3) = -1e10;
    end
else
    log_L_components(3) = -1e10;
    pred.cph = [];
end

% -------------------- 4) Group velocity --------------------
if isfield(data, 'ugr') && ~isempty(data.ugr.T) && pred_raw.ugr.status == 0
    pred.ugr = pred_raw.ugr.val(:);

    obs = data.ugr.obs(:);
    sig = data.ugr.sig(:);

    valid = isfinite(obs) & isfinite(sig) & (sig > 0) & isfinite(pred.ugr);
    n_valid = sum(valid);
    n_total = numel(obs);

    if n_valid > 0
        log_L_components(4) = gaussian_log_likelihood_hetero(obs(valid), pred.ugr(valid), sig(valid), sigma_e(4));
        if n_valid < n_total
            log_L_components(4) = log_L_components(4) * (n_valid / n_total);
        end
    else
        log_L_components(4) = -1e10;
    end
else
    log_L_components(4) = -1e10;
    pred.ugr = [];
end

log_L = sum(log_L_components);
if ~isfinite(log_L)
    log_L = -1e10;
end

end

% =======================================================================

function model = build_model_matrix(h, vs)

n_layers = 6;
model = zeros(n_layers, 4);

for i = 1:n_layers
    if i <= 5
        model(i, 1) = h(i);
    else
        model(i, 1) = 0;
    end
    model(i, 3) = vs(i);
    [model(i, 2), model(i, 4)] = brocher2005(vs(i));
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

n_hvsr = 0; n_ellip = 0; n_cph = 0; n_ugr = 0;
if isfield(data, 'hvsr') && isfield(data.hvsr, 'f')
    n_hvsr = sum(isfinite(data.hvsr.obs));
end
if isfield(data, 'ellip') && isfield(data.ellip, 'T')
    n_ellip = sum(isfinite(data.ellip.obs));
end
if isfield(data, 'cph') && isfield(data.cph, 'T')
    n_cph = sum(isfinite(data.cph.obs));
end
if isfield(data, 'ugr') && isfield(data.ugr, 'T')
    n_ugr = sum(isfinite(data.ugr.obs));
end

config.use_hvsr  = (n_hvsr  > 0);
config.use_ellip = (n_ellip > 0);
config.use_cph   = (n_cph   > 0);
config.use_ugr   = (n_ugr   > 0);

%% Parameters for HVf (your requested values)
config.hvsr.nf   = 100;
config.hvsr.fmin = 1/3;
config.hvsr.fmax = 25;
config.hvsr.nmr  = 10;
config.hvsr.nml  = 10;
config.hvsr.nks  = 1000; % keep your stability increase

config.ellip.nf   = 100;
config.disp.nmode = 0;

end

function fwd_data = build_forward_data(data)

fwd_data = struct();

if isfield(data, 'hvsr') && ~isempty(data.hvsr.f)
    fwd_data.hvsr.f   = data.hvsr.f(:);
    fwd_data.hvsr.obs = data.hvsr.obs(:);
    fwd_data.hvsr.sig = data.hvsr.sig(:);
end

if isfield(data, 'ellip') && ~isempty(data.ellip.T)
    fwd_data.ellip.T   = data.ellip.T(:);
    fwd_data.ellip.obs = data.ellip.obs(:);
    fwd_data.ellip.sig = data.ellip.sig(:);
end

if isfield(data, 'cph') && ~isempty(data.cph.T)
    fwd_data.cph.T   = data.cph.T(:);
    fwd_data.cph.obs = data.cph.obs(:);
    fwd_data.cph.sig = data.cph.sig(:);
end

if isfield(data, 'ugr') && ~isempty(data.ugr.T)
    fwd_data.ugr.T   = data.ugr.T(:);
    fwd_data.ugr.obs = data.ugr.obs(:);
    fwd_data.ugr.sig = data.ugr.sig(:);
end

end

function log_L = gaussian_log_likelihood_hetero(obs, pred, sig_obs, sigma_mult)
% Heteroscedastic Gaussian log-likelihood with sigma_mult scaling.
%
% Effective std per point: sigma_eff = sigma_mult * sig_obs

if isempty(pred) || any(~isfinite(pred)) || any(~isfinite(obs)) || any(~isfinite(sig_obs))
    log_L = -1e10;
    return;
end

sigma_mult = max(sigma_mult, 1e-6);
sig_obs = sig_obs(:);
obs     = obs(:);
pred    = pred(:);

sigma_eff = sigma_mult .* sig_obs;
if any(sigma_eff <= 0) || any(~isfinite(sigma_eff))
    log_L = -1e10;
    return;
end

r = (obs - pred) ./ sigma_eff;

% log N(obs | pred, sigma_eff^2)
log_L = -0.5 * numel(obs) * log(2*pi) ...
        - sum(log(sigma_eff)) ...
        - 0.5 * sum(r.^2);

end
