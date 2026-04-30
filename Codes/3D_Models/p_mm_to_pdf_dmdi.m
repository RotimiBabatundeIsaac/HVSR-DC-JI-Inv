function [pdf_mod] = p_mm_to_pdf_dmdi(vsta_mod, pdf_terp, min_mm_terp, dmm_di, nmm, nsta)
% Fast function to interpolate the pdf at each station based on chosen
% velocities (mm) at those stations (vsta_mod). 
% To avoid any sorting or excess math, I pick the two nearest velocity
% indecies of pdf_terp by converting vsta_mod directly into indecies. This
% is done using the ratio of how much v changes per index (dmm_di). 
%
% vsta_mod: (nsta x 1) Velocity at each station. Determines where to sample PDFs. 
% pdf_terp: (nmm x nsta) pdf at each station. Interpolated to a common mm. 
% min_mm_terp: (scalar) Smallest mm (velocity) used to make pdf_terp. 
% dmm_di: (scalar) How much does velocity (mm) change with each index in mm_terp? 
% nmm: (scalar) Number of mm or velocity in mm_terp. 
% nsta: (scalar) Number of stations. 
%
% brb2022.10.05

% For later interpolation, find two nearest indecies corresponding to a mm value. Use dmm_di rather than mm_terp. This should be much faster than using distances and sorting because we calculate indecies directly. 
best_ind = 1 + (vsta_mod - min_mm_terp)/dmm_di;   % Start at one. Find (non integer) closest index of mm_terp based on dindex/dmm
best_ind_flce = [floor(best_ind), ceil(best_ind)];   % Indecies before (floor) and after (ceil) the idealized non-integer index. 

% Handle bounds problems. 
best_ind_flce(best_ind_flce<1)=1;   % If lower than our bounds, just use our lowest part of pdf. That should be 0 anyway!
best_ind_flce(best_ind_flce>nmm) = nmm;   % If greater than our bounds, just use our last part of pdf. That should be 0 anyway!

% Convert nmm x nsta coordinates to linear coordinates, to do vectorized operations more easily. 
sta_hack = [[1:nsta];[1:nsta]]';   % Which station each model parameter corresponds, in best_ind_flce. % TODO can be taken outside function if really a little more speed. 
best_ind_linear = sub2ind([nmm, nsta], best_ind_flce, sta_hack);   % Linear (single dimension) coordinates to points of interest in pdf_terp

% Weights for interpolation
ind_weight = 1./([best_ind - best_ind_flce].^2);   % Convert difference between mm and the nearest two mm_terp points to a weight. 
ind_weight(or(isnan(ind_weight),ind_weight>99999)) = 99999;   % In case of 1./0. 
ind_weight = ind_weight ./ sum(ind_weight,2);   % Normalize, so multiply this to values at ind and this gives an average. 

% Do the interpolation
pdf_mod = sum(pdf_terp(best_ind_linear).*ind_weight,2);   % Interpolated pdf between two mm values

end