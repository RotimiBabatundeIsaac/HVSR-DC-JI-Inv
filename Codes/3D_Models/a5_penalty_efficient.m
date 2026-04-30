function [penalty, penalty_no_prior, roughness, pen_norm] = a5_penalty_efficient(vgrid,...
    pdf_terp, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay,...
    nearesti, weighti, min_mm_terp, dmm_di, nmm, nsta)
% Penalty function for surface inversion using fminunc
% Balances PDF likelihood at stations against spatial roughness
% Adapted from Brunsvik et al. for SESAME HBI

% Second-derivative roughness
dvdx2 = ((vgrid(1:end-2,:) - 2*vgrid(2:end-1,:) + vgrid(3:end,:))./dx2);
dvdy2 = ((vgrid(:,1:end-2) - 2*vgrid(:,2:end-1) + vgrid(:,3:end))./dy2);
roughness = sum(dvdx2 .^ 2, 'all') + sum(dvdy2 .^ 2, 'all');
roughness = roughness * rough_scale;

% Interpolate model values at station locations using precomputed weights
vsta = sum(vgrid(nearesti) .* weighti, 2);

% Look up PDF probability at each station's model value
pv_mod = p_mm_to_pdf_dmdi(vsta, pdf_terp, min_mm_terp, dmm_di, nmm, nsta);

% PDF penalty (negative sum = we want to maximize probability)
penalty_no_prior = -sum(pv_mod);

% Norm penalty for stability far from stations
dvgrid = vgrid - mean(vgrid, 'all');
dnorm = sqrt(dvgrid(:)' * dvgrid(:));
pen_norm = 1e-1 * dnorm;

penalty = penalty_no_prior + roughness + pen_norm;
end