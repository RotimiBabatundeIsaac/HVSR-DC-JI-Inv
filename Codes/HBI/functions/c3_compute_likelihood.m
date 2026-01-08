function [log_L, log_L_components, pred] = c3_compute_likelihood(theta, sigma_e, data, prior)
% Compute log-likelihood for all four data types
%
% Inputs:
%   theta   - model parameters [11 x 1]
%   sigma_e - noise standard deviations [4 x 1]
%   data    - structure with observed data for station
%   prior   - structure from c1_define_hbi_prior (contains paths)
%
% Outputs:
%   log_L            - total log-likelihood (scalar)
%   log_L_components - log-likelihood for each data type [4 x 1]
%   pred             - structure with predicted data

% Set up paths only once (use persistent variable)
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

% Extract parameters and build model matrix
h = theta(1:5);
vs = theta(6:11);
model = build_model_matrix(h, vs);

% Build config structure
config = build_forward_config(data);

% Build data structure for forward_model
fwd_data = build_forward_data(data);

% Run forward model
pred_raw = forward_model(model, fwd_data, config);

% Initialize outputs
log_L_components = zeros(4, 1);
pred = struct();

% 1. HVSR
if isfield(data, 'hvsr') && ~isempty(data.hvsr.f) && pred_raw.hvsr.status == 0
    pred.hvsr = pred_raw.hvsr.val;
    valid = isfinite(data.hvsr.obs) & isfinite(pred.hvsr);
    if sum(valid) > 0
        log_L_components(1) = gaussian_log_likelihood(data.hvsr.obs(valid), pred.hvsr(valid), sigma_e(1));
    else
        log_L_components(1) = -1e10;
    end
else
    log_L_components(1) = -1e10;
    pred.hvsr = [];
end

% 2. Ellipticity
if isfield(data, 'ellip') && ~isempty(data.ellip.T) && pred_raw.ellip.status == 0
    pred.ellip = pred_raw.ellip.val;
    valid = isfinite(data.ellip.obs) & isfinite(pred.ellip);
    if sum(valid) > 0
        log_L_components(2) = gaussian_log_likelihood(data.ellip.obs(valid), pred.ellip(valid), sigma_e(2));
    else
        log_L_components(2) = -1e10;
    end
else
    log_L_components(2) = -1e10;
    pred.ellip = [];
end

% 3. Phase velocity
if isfield(data, 'cph') && ~isempty(data.cph.T) && pred_raw.cph.status == 0
    pred.cph = pred_raw.cph.val;
    valid = isfinite(data.cph.obs) & isfinite(pred.cph);
    if sum(valid) > 0
        log_L_components(3) = gaussian_log_likelihood(data.cph.obs(valid), pred.cph(valid), sigma_e(3));
    else
        log_L_components(3) = -1e10;
    end
else
    log_L_components(3) = -1e10;
    pred.cph = [];
end

% 4. Group velocity
if isfield(data, 'ugr') && ~isempty(data.ugr.T) && pred_raw.ugr.status == 0
    pred.ugr = pred_raw.ugr.val;
    valid = isfinite(data.ugr.obs) & isfinite(pred.ugr);
    if sum(valid) > 0
        log_L_components(4) = gaussian_log_likelihood(data.ugr.obs(valid), pred.ugr(valid), sigma_e(4));
    else
        log_L_components(4) = -1e10;
    end
else
    log_L_components(4) = -1e10;
    pred.ugr = [];
end

% Total log-likelihood
log_L = sum(log_L_components);

% Handle non-finite values
if ~isfinite(log_L)
    log_L = -1e10;
end

end


function model = build_model_matrix(h, vs)
% Build model matrix [h, Vp, Vs, rho] from parameters
%
% Inputs:
%   h  - layer thicknesses [5 x 1] (km)
%   vs - shear velocities [6 x 1] (km/s)
%
% Output:
%   model - [6 x 4] matrix: [thickness, Vp, Vs, density]

n_layers = 6;
model = zeros(n_layers, 4);

for i = 1:n_layers
    if i <= 5
        model(i, 1) = h(i);
    else
        model(i, 1) = 0;  % half-space
    end
    model(i, 3) = vs(i);
    [model(i, 2), model(i, 4)] = brocher2005(vs(i));
end

end


function [vp, rho] = brocher2005(vs)
% Empirical relations from Brocher (2005)
%
% Input:
%   vs - shear wave velocity (km/s)
%
% Output:
%   vp  - P-wave velocity (km/s)
%   rho - density (g/cm^3)

vp = 0.9409 + 2.0947*vs - 0.8206*vs^2 + 0.2683*vs^3 - 0.0251*vs^4;
rho = 1.6612*vp - 0.4721*vp^2 + 0.0671*vp^3 - 0.0043*vp^4 + 0.000106*vp^5;

% Ensure physical bounds
vp = max(vp, vs * 1.5);
rho = max(rho, 1.5);

end


function config = build_forward_config(data)
% Build config structure for forward_model

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

config.use_hvsr = (n_hvsr > 0);
config.use_ellip = (n_ellip > 0);
config.use_cph = (n_cph > 0);
config.use_ugr = (n_ugr > 0);

config.hvsr.nf = 200;
config.hvsr.fmin = 0.1;
config.hvsr.fmax = 30;
config.hvsr.nmr = 10;
config.hvsr.nml = 10;
config.hvsr.nks = 500;
config.ellip.nf = 100;
config.disp.nmode = 0;

end


function fwd_data = build_forward_data(data)
% Build data structure for forward_model

fwd_data = struct();

if isfield(data, 'hvsr') && ~isempty(data.hvsr.f)
    fwd_data.hvsr.f = data.hvsr.f(:);
    fwd_data.hvsr.obs = data.hvsr.obs(:);
    fwd_data.hvsr.sig = data.hvsr.sig(:);
end

if isfield(data, 'ellip') && ~isempty(data.ellip.T)
    fwd_data.ellip.T = data.ellip.T(:);
    fwd_data.ellip.obs = data.ellip.obs(:);
    fwd_data.ellip.sig = data.ellip.sig(:);
end

if isfield(data, 'cph') && ~isempty(data.cph.T)
    fwd_data.cph.T = data.cph.T(:);
    fwd_data.cph.obs = data.cph.obs(:);
    fwd_data.cph.sig = data.cph.sig(:);
end

if isfield(data, 'ugr') && ~isempty(data.ugr.T)
    fwd_data.ugr.T = data.ugr.T(:);
    fwd_data.ugr.obs = data.ugr.obs(:);
    fwd_data.ugr.sig = data.ugr.sig(:);
end

end


function log_L = gaussian_log_likelihood(obs, pred, sigma)
% Compute Gaussian log-likelihood
%
% Inputs:
%   obs   - observed data [N x 1]
%   pred  - predicted data [N x 1]
%   sigma - noise standard deviation (scalar)
%
% Output:
%   log_L - log-likelihood (scalar)

if isempty(pred) || any(isnan(pred)) || any(~isfinite(pred))
    log_L = -1e10;
    return;
end

N = length(obs);
residuals = obs(:) - pred(:);
chi2 = sum(residuals.^2);

log_L = -0.5 * N * log(2 * pi) - N * log(sigma) - 0.5 * chi2 / (sigma^2);

end