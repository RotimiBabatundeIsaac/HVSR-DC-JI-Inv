function Sigma_theta_new = c7_gibbs_update_Sigma_theta(theta, mu_theta, prior)
% Gibbs sample for Sigma_theta from conditional Inverse-Wishart distribution
%
% For single observation (one station):
%   Sigma_theta | theta, mu_theta ~ IW(nu_post, Psi_post)
%   where:
%     nu_post = nu_0 + 1
%     Psi_post = Psi_0 + (theta - mu_theta)(theta - mu_theta)'
%
% Inputs:
%   theta    - current model parameters [11 x 1]
%   mu_theta - current mean vector [11 x 1]
%   prior    - structure from c1_define_hbi_prior
%
% Output:
%   Sigma_theta_new - sampled covariance matrix [11 x 11]

n = length(theta);

% Posterior parameters
nu_post = prior.nu_0 + 1;

diff = theta - mu_theta;
Psi_post = prior.Psi_0 + (diff * diff');

% Sample from Inverse-Wishart(nu_post, Psi_post)
Sigma_theta_new = sample_inverse_wishart(nu_post, Psi_post);

% Ensure positive definiteness
[~, p] = chol(Sigma_theta_new);
if p > 0
    warning('Sampled Sigma_theta not positive definite, using regularization');
    Sigma_theta_new = Sigma_theta_new + 1e-6 * eye(n);
end

end


function S = sample_inverse_wishart(nu, Psi)
% Sample from Inverse-Wishart distribution IW(nu, Psi)
%
% Method: If X ~ Wishart(nu, Psi^-1), then X^-1 ~ IW(nu, Psi)
%
% Inputs:
%   nu  - degrees of freedom (scalar)
%   Psi - scale matrix [n x n]
%
% Output:
%   S - sampled covariance matrix [n x n]

n = size(Psi, 1);

% Cholesky of Psi^-1
[L_Psi, p] = chol(Psi, 'lower');
if p > 0
    % Psi not positive definite, add regularization
    Psi = Psi + 1e-6 * eye(n);
    L_Psi = chol(Psi, 'lower');
end
L_Psi_inv = L_Psi \ eye(n);

% Sample from Wishart using Bartlett decomposition
A = zeros(n, n);
for i = 1:n
    A(i, i) = sqrt(chi2rnd(nu - i + 1));
    for j = 1:i-1
        A(i, j) = randn();
    end
end

% Wishart sample: W = L * A * A' * L'
L_W = L_Psi_inv * A;
W = L_W * L_W';

% Inverse-Wishart sample
S = W \ eye(n);

% Ensure symmetry
S = (S + S') / 2;

end