function prior = c1_define_hbi_prior()
% Define all hyperprior parameters, bounds, and constraints for HBI

prior.n_params = 11;
prior.use_hierarchical = true;

prior.param_names = {'h1','h2','h3','h4','h5','Vs1','Vs2','Vs3','Vs4','Vs5','Vs6'};
prior.param_units = {'km','km','km','km','km','km/s','km/s','km/s','km/s','km/s','km/s'};

prior.idx.h = 1:5;
prior.idx.vs = 6:11;

prior.enforce_monotonicity = true;
prior.monotonic_start_idx = 6;

% -------------------- BOUNDS --------------------
% Keep deep structure flexible for dispersion, but tighten shallow space
% to make the HVSR/ellip peak region easier to find.

% Thickness bounds (km)
h_min = [0.005, 0.02, 0.08,  3.0, 10.0]';   % keep h4/h5 like your current
h_max = [0.080, 0.30, 1.50, 25.0, 45.0]';   % tighten h1-h3

% Vs bounds (km/s)
% Tighten shallow Vs to encourage strong contrasts that control HVSR peak shape
vs_min = [0.12, 0.25, 0.60, 2.50, 2.80, 3.20]';
vs_max = [0.90, 1.60, 2.60, 4.20, 4.60, 5.20]';

prior.hard_min = [h_min; vs_min];
prior.hard_max = [h_max; vs_max];

% -------------------- INITIAL MODEL --------------------
% Start closer to a sharp ~8–10 Hz HVSR peak mechanism: thin/slow cap over faster layer.
prior.mu_theta_init = [ ...
    0.020, 0.120, 0.600, 10.00, 24.00, ...  % h1-h5 (km)
    0.35,  0.75,  1.40,  3.53,  3.73, 4.40  ... % Vs1-Vs6 (km/s)
]';

% -------------------- HIERARCHICAL PRIORS --------------------
param_range = prior.hard_max - prior.hard_min;
prior.prior_std = 0.50 * param_range;

prior.nu_0 = 13;
prior.Psi_0 = diag(prior.prior_std.^2);

% -------------------- NOISE PRIORS --------------------
% Prevent sigma_hvsr/sigma_ellip from inflating too easily.
prior.noise.alpha_0 = 2.0;
prior.noise.beta_0  = 0.01;

prior.noise.names = {'sigma_hvsr','sigma_ellip','sigma_cph','sigma_ugr'};
prior.noise.n_noise = 4;

prior.noise.min  = [0.01, 0.01, 0.01, 0.01]';
prior.noise.max  = [0.80, 0.60, 0.50, 0.50]';   % tightened HVSR/ellip caps
prior.noise.init = [0.15, 0.08, 0.05, 0.08]';   % start with stronger HVSR/ellip weighting

% -------------------- PROPOSAL STEP SIZES --------------------
% Much smaller steps for h1-h3 and Vs1-Vs3 to improve acceptance and exploration near HVSR peak.
prior.step_size = [ ...
    0.003, 0.010, 0.030, 1.0,  1.5, ...   % h1-h5
    0.03,  0.05,  0.08,  0.08, 0.08, 0.06 ...  % Vs1-Vs6
]';

% -------------------- MCMC SETTINGS --------------------
prior.mcmc.n_iterations = 10000;
prior.mcmc.burn_in = 1000;
prior.mcmc.thin = 50;   % 200 is very aggressive; reduces posterior resolution


% Paths to forward modeling software (leave unchanged)
prior.paths.base = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion';
prior.paths.cps_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_v3.30');
prior.paths.hv_dfa_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_HV_DFA');

end
