function [theta_prop, is_valid, perturb_info] = c5_perturb_theta(theta, prior)
%c5_perturb_theta Perturb model parameters for Metropolis-Hastings proposal

n_params = length(theta);

step_sizes = get_step_sizes(prior, n_params);

perturbations = randn(n_params, 1) .* step_sizes;
theta_prop = theta + perturbations;

[lb, ub] = get_bounds(prior);

is_valid = true;
if any(theta_prop < lb) || any(theta_prop > ub)
    is_valid = false;
end

if is_valid
    if isfield(prior, 'enforce_monotonicity') && ~prior.enforce_monotonicity
    else
        vs_prop = get_vs(theta_prop, prior);
        if any(diff(vs_prop) < 0)
            is_valid = false;
        end
    end
end

perturb_info.step_sizes = step_sizes;
perturb_info.perturbations = perturbations;

end

function step_sizes = get_step_sizes(prior, n_params)
    if isfield(prior, 'mcmc') && isfield(prior.mcmc, 'step_size')
        step_sizes = prior.mcmc.step_size(:);
    elseif isfield(prior, 'step_size')
        step_sizes = prior.step_size(:);
    else
        n_h = 5;
        n_vs = 6;
        step_sizes = [0.5*ones(n_h,1); 0.1*ones(n_vs,1)];
    end
    
    if length(step_sizes) ~= n_params
        step_sizes = ones(n_params, 1) * 0.1;
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

function vs = get_vs(theta, prior)
    if isfield(prior, 'idx') && isfield(prior.idx, 'vs')
        vs = theta(prior.idx.vs);
    elseif isfield(prior, 'n_layers')
        n_h = prior.n_layers - 1;
        vs = theta(n_h+1:end);
    else
        n_params = length(theta);
        n_vs = 6;
        n_h = n_params - n_vs;
        vs = theta(n_h+1:end);
    end
end