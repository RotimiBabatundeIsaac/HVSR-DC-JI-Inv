function [theta_new, is_valid, idx_perturbed] = c5_perturb_theta(theta, prior)
% Propose MH perturbation with constraint checking
%
% Inputs:
%   theta - current model parameters [11 x 1]
%   prior - structure from c1_define_hbi_prior
%
% Outputs:
%   theta_new     - proposed model parameters [11 x 1]
%   is_valid      - boolean, true if proposal satisfies all constraints
%   idx_perturbed - index of the perturbed parameter

n_params = prior.n_params;

% Randomly select one parameter to perturb
idx_perturbed = randi(n_params);

% Copy current theta
theta_new = theta;

% Propose perturbation using Gaussian random walk
step = prior.step_size(idx_perturbed) * randn();
theta_new(idx_perturbed) = theta(idx_perturbed) + step;

% Check constraints
is_valid = check_constraints(theta_new, prior);

end


function is_valid = check_constraints(theta, prior)
% Check all constraints on theta
%
% Constraints:
%   1. Within hard bounds
%   2. Velocity monotonicity: Vs1 <= Vs2 <= Vs3 <= Vs4 <= Vs5 <= Vs6
%   3. Positivity of thicknesses

is_valid = true;

% Check hard bounds
if any(theta < prior.hard_min) || any(theta > prior.hard_max)
    is_valid = false;
    return;
end

% Check velocity monotonicity
vs = theta(6:11);
if any(diff(vs) < 0)
    is_valid = false;
    return;
end

% Check positivity of thicknesses
h = theta(1:5);
if any(h <= 0)
    is_valid = false;
    return;
end

end