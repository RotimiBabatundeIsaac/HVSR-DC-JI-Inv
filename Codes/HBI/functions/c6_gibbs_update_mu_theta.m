function mu_theta_new = c6_gibbs_update_mu_theta(theta, Sigma_theta, prior)
% Gibbs sample for mu_theta from conditional Normal distribution
%
% For single observation (one station):
%   mu_theta | theta, Sigma_theta ~ N(theta, Sigma_theta)
%
% Inputs:
%   theta       - current model parameters [11 x 1]
%   Sigma_theta - current covariance matrix [11 x 11]
%   prior       - structure from c1_define_hbi_prior
%
% Output:
%   mu_theta_new - sampled mean vector [11 x 1]

n = length(theta);

% Cholesky decomposition for sampling
[L, p] = chol(Sigma_theta, 'lower');

if p > 0
    % Sigma_theta not positive definite, return theta as fallback
    warning('Sigma_theta not positive definite in mu_theta update');
    mu_theta_new = theta;
    return;
end

% Sample from N(theta, Sigma_theta)
z = randn(n, 1);
mu_theta_new = theta + L * z;

% Enforce hard bounds on mu_theta
mu_theta_new = max(mu_theta_new, prior.hard_min);
mu_theta_new = min(mu_theta_new, prior.hard_max);

end