function [log_prior, is_valid] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior)
%c4_compute_log_prior_theta Compute log-prior for model parameters theta

n_params = length(theta);

[lb, ub] = get_bounds(prior);
is_valid = check_constraints(theta, lb, ub, prior);

if ~is_valid
    log_prior = -Inf;
    return;
end

use_hierarchical = get_hierarchical_flag(prior);

if use_hierarchical
    diff = theta - mu_theta;
    log_prior = -0.5 * n_params * log(2*pi) ...
                - 0.5 * log(det(Sigma_theta)) ...
                - 0.5 * (diff' / Sigma_theta * diff);
else
    log_prior = 0;
    for i = 1:n_params
        log_prior = log_prior - log(ub(i) - lb(i));
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

function use_hier = get_hierarchical_flag(prior)
    if isfield(prior, 'use_hierarchical_prior')
        use_hier = prior.use_hierarchical_prior;
    elseif isfield(prior, 'use_hierarchical')
        use_hier = prior.use_hierarchical;
    else
        use_hier = true;
    end
end

function is_valid = check_constraints(theta, lb, ub, prior)
    is_valid = true;
    
    if any(theta < lb) || any(theta > ub)
        is_valid = false;
        return;
    end
    
    if isfield(prior, 'enforce_monotonicity') && ~prior.enforce_monotonicity
        return;
    end
    
    vs = get_vs(theta, prior);
    if any(diff(vs) < 0)
        is_valid = false;
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