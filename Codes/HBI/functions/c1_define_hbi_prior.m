function prior = c1_define_hbi_prior()
% Define all hyperprior parameters, bounds, and constraints for HBI
%
% Output:
%   prior - structure containing all prior specifications

% Number of model parameters
prior.n_params = 11;



%we need control to switch between guassian and non guassian 


prior.use_hierarchical =true; %false;  % true = Gaussian HBI, false = uniform prior MCMC


% Parameter names for reference
prior.param_names = {'h1', 'h2', 'h3', 'h4', 'h5', ...
                     'Vs1', 'Vs2', 'Vs3', 'Vs4', 'Vs5', 'Vs6'};

% Parameter units
prior.param_units = {'km', 'km', 'km', 'km', 'km', ...
                     'km/s', 'km/s', 'km/s', 'km/s', 'km/s', 'km/s'};

% Hard bounds [min, max] for each parameter
%   h1    h2    h3    h4    h5   Vs1   Vs2   Vs3   Vs4   Vs5   Vs6
prior.hard_min = [0.005, 0.02, 0.10,  3.0, 10.0, 0.10, 0.20, 0.40, 2.50, 2.80, 3.20]';
prior.hard_max = [0.300, 0.80, 3.00, 25.0, 45.0, 1.80, 2.20, 3.00, 4.20, 4.60, 5.20]';

% Initial mean values for mu_theta
% Based on test model that fits W21 data well
prior.mu_theta_init = [0.03, 0.17, 0.53, 10.27, 24.0, ...
                       0.50, 0.80, 0.99, 3.53, 3.73, 4.40]';

% Prior standard deviations (used to construct Psi_0)
% Set to 50% of the parameter range following Baise et al.
param_range = prior.hard_max - prior.hard_min;
prior.prior_std = 0.75 * param_range;

% Inverse-Wishart hyperprior for Sigma_theta
prior.nu_0 = 13;  % Degrees of freedom (weakly informative, nu_0 > n_params + 1)
prior.Psi_0 = diag(prior.prior_std.^2);  % Scale matrix (diagonal)

% Inverse-Gamma hyperprior for noise variances
prior.noise.alpha_0 = 1.0;  % Shape parameter
prior.noise.beta_0 = 0.01;  % Scale parameter

% Noise parameter specifications
prior.noise.names = {'sigma_hvsr', 'sigma_ellip', 'sigma_cph', 'sigma_ugr'};
prior.noise.n_noise = 4;
prior.noise.min = [0.01, 0.01, 0.01, 0.01]';
prior.noise.max = [2.0, 2.0, 0.5, 0.5]';
prior.noise.init = [0.30, 0.10, 0.05, 0.08]';

% MH step sizes for theta perturbation (tuned for ~30% acceptance)
prior.step_size = [0.01, 0.03, 0.10, 1.5, 2.0, ...
                   0.08, 0.10, 0.12, 0.10, 0.10, 0.08]';

% MCMC settings
prior.mcmc.n_iterations = 1000;%number of mcmc iterations 
prior.mcmc.burn_in = 200;% Discard the n number of starting burn
prior.mcmc.thin = 2; % how each should be saved to the posterior 

% Paths to forward modeling software
prior.paths.base = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion';
prior.paths.cps_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_v3.30');
prior.paths.hv_dfa_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_HV_DFA');

end