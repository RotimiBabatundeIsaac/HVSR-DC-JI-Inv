function sigma_e_new = c8_gibbs_update_sigma_e(pred, data, prior)
% Gibbs sample for noise parameters from conditional Inverse-Gamma distribution
%
% For each data type d with N_d observations:
%   sigma_d^2 | data, theta ~ IG(alpha_post, beta_post)
%   where:
%     chi2_d = sum((obs_d - pred_d)^2)
%     alpha_post = alpha_0 + N_d/2
%     beta_post = beta_0 + chi2_d/2
%
% Inputs:
%   pred  - structure with predicted data from forward model
%   data  - structure with observed data
%   prior - structure from c1_define_hbi_prior
%
% Output:
%   sigma_e_new - sampled noise standard deviations [4 x 1]

sigma_e_new = prior.noise.init;  % Initialize with prior values

alpha_0 = prior.noise.alpha_0;
beta_0 = prior.noise.beta_0;

% 1. HVSR
if isfield(data, 'hvsr') && ~isempty(data.hvsr.f) && ~isempty(pred.hvsr)
    valid = isfinite(data.hvsr.obs) & isfinite(pred.hvsr);
    if sum(valid) > 0
        residuals = data.hvsr.obs(valid) - pred.hvsr(valid);
        chi2 = sum(residuals.^2);
        N_d = sum(valid);
        alpha_post = alpha_0 + N_d / 2;
        beta_post = beta_0 + chi2 / 2;
        sigma2 = sample_inverse_gamma(alpha_post, beta_post);
        sigma_e_new(1) = sqrt(sigma2);
    end
end

% 2. Ellipticity
if isfield(data, 'ellip') && ~isempty(data.ellip.T) && ~isempty(pred.ellip)
    valid = isfinite(data.ellip.obs) & isfinite(pred.ellip);
    if sum(valid) > 0
        residuals = data.ellip.obs(valid) - pred.ellip(valid);
        chi2 = sum(residuals.^2);
        N_d = sum(valid);
        alpha_post = alpha_0 + N_d / 2;
        beta_post = beta_0 + chi2 / 2;
        sigma2 = sample_inverse_gamma(alpha_post, beta_post);
        sigma_e_new(2) = sqrt(sigma2);
    end
end

% 3. Phase velocity
if isfield(data, 'cph') && ~isempty(data.cph.T) && ~isempty(pred.cph)
    valid = isfinite(data.cph.obs) & isfinite(pred.cph);
    if sum(valid) > 0
        residuals = data.cph.obs(valid) - pred.cph(valid);
        chi2 = sum(residuals.^2);
        N_d = sum(valid);
        alpha_post = alpha_0 + N_d / 2;
        beta_post = beta_0 + chi2 / 2;
        sigma2 = sample_inverse_gamma(alpha_post, beta_post);
        sigma_e_new(3) = sqrt(sigma2);
    end
end

% 4. Group velocity
if isfield(data, 'ugr') && ~isempty(data.ugr.T) && ~isempty(pred.ugr)
    valid = isfinite(data.ugr.obs) & isfinite(pred.ugr);
    if sum(valid) > 0
        residuals = data.ugr.obs(valid) - pred.ugr(valid);
        chi2 = sum(residuals.^2);
        N_d = sum(valid);
        alpha_post = alpha_0 + N_d / 2;
        beta_post = beta_0 + chi2 / 2;
        sigma2 = sample_inverse_gamma(alpha_post, beta_post);
        sigma_e_new(4) = sqrt(sigma2);
    end
end

% Enforce bounds
sigma_e_new = max(sigma_e_new, prior.noise.min);
sigma_e_new = min(sigma_e_new, prior.noise.max);

end


function x = sample_inverse_gamma(alpha, beta)
% Sample from Inverse-Gamma distribution IG(alpha, beta)
%
% Method: If Y ~ Gamma(alpha, 1/beta), then 1/Y ~ IG(alpha, beta)
%
% Inputs:
%   alpha - shape parameter
%   beta  - scale parameter
%
% Output:
%   x - sampled value

y = gamrnd(alpha, 1/beta);
x = 1 / y;

end