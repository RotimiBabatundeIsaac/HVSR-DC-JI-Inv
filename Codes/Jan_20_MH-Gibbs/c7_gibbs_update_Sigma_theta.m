function Sigma_theta_new = c7_gibbs_update_Sigma_theta(theta, mu_theta, prior)
%c7_gibbs_update_Sigma_theta Gibbs update for hyperparameter Sigma_theta

n_params = length(theta);

if isfield(prior, 'hyper') && isfield(prior.hyper, 'nu0')
    nu0 = prior.hyper.nu0;
    Psi0 = prior.hyper.Psi0;
elseif isfield(prior, 'nu_0')
    nu0 = prior.nu_0;
    Psi0 = prior.Psi_0;
else
    [lb, ub] = get_bounds(prior);
    nu0 = n_params + 2;
    scale = (ub - lb) / 6;
    Psi0 = diag(scale.^2) * nu0;
end

diff = theta - mu_theta;
S = diff * diff';

nu_post = nu0 + 1;
Psi_post = Psi0 + S;

Sigma_theta_new = iwishrnd(Psi_post, nu_post);

min_var = 1e-6;
for i = 1:n_params
    if Sigma_theta_new(i,i) < min_var
        Sigma_theta_new(i,i) = min_var;
    end
end

end

function [lb, ub] = get_bounds(prior)
    if isfield(prior, 'bounds')
        lb = prior.bounds.lower;
        ub = prior.bounds.upper;
    elseif isfield(prior, 'hard_min')
        lb = prior.hard_min;
        ub = prior.hard_max;
    else
        error('No bounds found in prior structure');
    end
end