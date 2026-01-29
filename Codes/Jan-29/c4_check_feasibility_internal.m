function [is_feasible, violation] = c4_check_feasibility_internal(theta, prior)

is_feasible = true;
violation = '';

if isempty(theta)
    is_feasible = false;
    violation = 'theta is empty';
    return;
end

if any(~isfinite(theta))
    is_feasible = false;
    idx = find(~isfinite(theta), 1);
    violation = sprintf('Parameter %d is not finite: %.4f', idx, theta(idx));
    return;
end

lb = prior.bounds.lower;
ub = prior.bounds.upper;

for i = 1:length(theta)
    if theta(i) < lb(i)
        is_feasible = false;
        violation = sprintf('Parameter %d below lower bound: %.6f < %.6f', i, theta(i), lb(i));
        return;
    end
    if theta(i) > ub(i)
        is_feasible = false;
        violation = sprintf('Parameter %d above upper bound: %.6f > %.6f', i, theta(i), ub(i));
        return;
    end
end

h = theta(prior.idx.h);
for i = 1:length(h)
    if h(i) <= 0
        is_feasible = false;
        violation = sprintf('Thickness h%d non-positive: %.6f', i, h(i));
        return;
    end
end

if prior.enforce_monotonicity
    vs = theta(prior.idx.vs);
    for i = 1:(length(vs)-1)
        if vs(i+1) < vs(i)
            is_feasible = false;
            violation = sprintf('Vs monotonicity violated: Vs%d=%.4f > Vs%d=%.4f', i, vs(i), i+1, vs(i+1));
            return;
        end
    end
end

end