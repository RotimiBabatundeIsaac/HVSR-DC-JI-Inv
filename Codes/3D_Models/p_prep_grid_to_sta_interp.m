function [fhand_vec, fhand_mat, grid_terp, nearesti, weighti...
    ] = p_prep_grid_to_sta_interp(xgrid, ygrid, vgrid, stax, stay)
%% Efficient interpolation from grid to stations. 
% Returns two function handles. Should make things easier than keeping
% track of a bunch of different weights and stuff. 
% 
% fhand_vec is vectorized operation. Provide it vgrid. It gives
% interpolated velocities. 
%
% fhand_mat is matrix-math operation. Give it (flattened) vgrid_r. Also
% given interpolated velocities. 
% 
% I suspect the matrix version is faster for smaller data/problems, and
% vector version faster for larger memory problems. 
% brb2022.10.05
nsta = length(stax); 

% Determine which nodes are closest to each station and find their weights.
ngrid = (size(xgrid,1)*size(xgrid,2)); 
vgrid_r = vgrid(:); % vgrid "reshaped" to 1d

nearesti = zeros(nsta, 4); % Four corners linear indices. 
weighti  = zeros(nsta, 4); % Weights to the four corners for interpolation. 
grid_terp = zeros(ngrid, nsta); % Matrix math version of interpolation. 

for ista = 1:nsta; 
    staxi = stax(ista); 
    stayi = stay(ista); 

    to_east  = xgrid >= staxi; 
    to_west  = xgrid <  staxi; 
    to_north = ygrid >= stayi; 
    to_south = ygrid <  stayi;

    xdist = xgrid - staxi;  
    ydist = ygrid - stayi; 
    tdist = sqrt(xdist.^2 + ydist.^2); 

    [easti,  ~] = find(to_east ); 
    [westi,  ~] = find(to_west ); 
    [~, northi] = find(to_north); 
    [~, southi] = find(to_south); 
    easti = min(easti); 
    westi = max(westi); 
    northi = min(northi); 
    southi = max(southi); 

    nei = sub2ind([size(vgrid,1), size(vgrid,2)], easti, northi); 
    nwi = sub2ind([size(vgrid,1), size(vgrid,2)], westi, northi); 
    sei = sub2ind([size(vgrid,1), size(vgrid,2)], easti, southi); 
    swi = sub2ind([size(vgrid,1), size(vgrid,2)], westi, southi); 

    box_map = [nei, nwi, sei, swi]; 

    nearesti(ista, :) = box_map; 

    dist_to_crnr = tdist(box_map); 
    weight_to_crnr = 1./dist_to_crnr; % TODO think about if this is the best interpolation method for this. 
    weight_to_crnr = weight_to_crnr ./ sum(weight_to_crnr); 
    
    if any(weight_to_crnr) == inf; 
        error('Should handle station being right on a node'); 
    end

    weighti(ista,:) = weight_to_crnr; 
    grid_terp(box_map,ista) = weight_to_crnr; 

end

if any(isnan(grid_terp)); 
    error('There are nans in grid_terp. Figure out a solution')
end

fhand_mat = @(vgrid_r)(vgrid_r' * grid_terp)';
fhand_vec = @(vgrid)sum(vgrid(nearesti) .* weighti,2);

% The rest is for testing. 
% % % 
% % % % Now can matrix multiply vgrid_r' * weight_to_crnr to interpolate velocity at
% % % % each station. Do not need to worry about distances again. 
% % % % I thnk also gridi_r needs to be passed to flatten vgrid
% % % vsta_terp_mat = (vgrid_r' * grid_terp)'; % Matrix multiplication of getting velocities. 
% % % vsta_terp_vec = sum(vgrid(nearesti) .* weighti,2); % Vectorized way of getting velocities. Less arithmatic. 
% % % % fhand_matmult = @()(vgrid_r' * grid_terp)';
% % % % fhand_matmult = @()(vgrid(gridi_r)' * grid_terp)';
% % % 
% % % %%% Now I just pass vgrid_r as argument and grid_terp as always present
% % % %%% argument
% % % 
% % % fhand_mat_timeit = @()(vgrid_r' * grid_terp)'; 
% % % fhand_vec_timeit = @()sum(vgrid(nearesti) .* weighti,2);
% % % timeit(fhand_mat_timeit) % fhand_mat might be a little faster, but if grid becomes extremely large or stations become very large, the number of operations grows by n or m squared
% % % timeit(fhand_vec_timeit) % Tiny bit slower with few stations or grid points. Maybe will be more efficient with more data. 
end