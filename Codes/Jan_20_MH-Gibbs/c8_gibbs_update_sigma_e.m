function sigma_e_new = c8_gibbs_update_sigma_e(pred, data, prior)
% Gibbs sample noise std parameters sigma_e for each datatype using IG update.
% Uses normalized residuals if data.(field).sig exists.

sigma_e_new = prior.noise.init(:);

alpha0 = prior.noise.alpha_0;
beta0  = prior.noise.beta_0;

alpha0 = expand_to_4(alpha0);
beta0  = expand_to_4(beta0);

sigma_e_new(1) = update_one_sigma(alpha0(1), beta0(1), ...
    get_obs(data,'hvsr'), get_pred(pred,'hvsr'), get_sig(data,'hvsr'));

sigma_e_new(2) = update_one_sigma(alpha0(2), beta0(2), ...
    get_obs(data,'ellip'), get_pred(pred,'ellip'), get_sig(data,'ellip'));

sigma_e_new(3) = update_one_sigma(alpha0(3), beta0(3), ...
    get_obs(data,'cph'), get_pred(pred,'cph'), get_sig(data,'cph'));

sigma_e_new(4) = update_one_sigma(alpha0(4), beta0(4), ...
    get_obs(data,'ugr'), get_pred(pred,'ugr'), get_sig(data,'ugr'));

% Enforce bounds
sigma_e_new = max(sigma_e_new, prior.noise.min(:));
sigma_e_new = min(sigma_e_new, prior.noise.max(:));

end

% ---------------- helpers ----------------

function s = update_one_sigma(alpha0, beta0, obs, prd, sig)
% Return a sampled sigma (std). If no valid data, return NaN and caller keeps init.
s = NaN;

if isempty(obs) || isempty(prd)
    return;
end

obs = obs(:); prd = prd(:);

% Ensure same length (guard against forward-model grid differences)
n = min(numel(obs), numel(prd));
obs = obs(1:n);
prd = prd(1:n);

valid = isfinite(obs) & isfinite(prd);

if ~isempty(sig)
    sig = sig(:);
    sig = sig(1:min(numel(sig), n));
    valid = valid & isfinite(sig) & (sig > 0);

    if ~any(valid); return; end

    % prevent tiny sig from exploding weights
    smin = prctile(sig(valid), 5);
    sig_eff = max(sig(valid), max(1e-6, 0.2*smin));

    r = (obs(valid) - prd(valid)) ./ sig_eff;
else
    if ~any(valid); return; end
    r = (obs(valid) - prd(valid));
end

chi2 = sum(r.^2);
Nd   = sum(valid);

alpha_post = alpha0 + Nd/2;
beta_post  = beta0  + chi2/2;

sigma2 = sample_inverse_gamma(alpha_post, beta_post);
s = sqrt(sigma2);
end

function obs = get_obs(data, field)
obs = [];
if isfield(data, field) && isfield(data.(field),'obs')
    obs = data.(field).obs;
end
end

function prd = get_pred(pred, field)
prd = [];
if isfield(pred, field)
    prd = pred.(field);
end
end

function sig = get_sig(data, field)
sig = [];
if isfield(data, field) && isfield(data.(field),'sig')
    sig = data.(field).sig;
end
end

function x = sample_inverse_gamma(alpha, beta)
% If Y ~ Gamma(alpha, 1/beta), then 1/Y ~ IG(alpha, beta)
y = gamrnd(alpha, 1/beta);
x = 1 ./ y;
end

function v = expand_to_4(x)
if isscalar(x)
    v = repmat(x, 4, 1);
else
    v = x(:);
    if numel(v) ~= 4
        error('prior.noise.alpha_0 and beta_0 must be scalar or length 4.');
    end
end
end
