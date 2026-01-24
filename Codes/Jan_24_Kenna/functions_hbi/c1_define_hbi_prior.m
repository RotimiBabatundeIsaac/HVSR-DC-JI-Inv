function prior = c1_define_hbi_prior()

prior.n_params = 11;
prior.param_names = {'h1','h2','h3','h4','h5','Vs1','Vs2','Vs3','Vs4','Vs5','Vs6'};
prior.param_units = {'km','km','km','km','km','km/s','km/s','km/s','km/s','km/s','km/s'};
prior.idx.h = 1:5;
prior.idx.vs = 6:11;

prior.enforce_monotonicity = true;
prior.monotonic_start_idx = 6;

h_min = [0.005; 0.02; 0.08; 3.0; 10.0];
h_max = [0.080; 0.30; 1.50; 25.0; 45.0];

vs_min = [0.12; 0.25; 0.60; 2.50; 2.80; 3.20];
vs_max = [0.90; 1.60; 2.60; 4.20; 4.60; 5.20];

prior.bounds.lower = [h_min; vs_min];
prior.bounds.upper = [h_max; vs_max];

prior.noise.n_noise = 4;
prior.noise.names = {'sigma_hvsr','sigma_ellip','sigma_cph','sigma_ugr'};
prior.noise.alpha_0 = [2.0; 2.0; 2.0; 2.0];
prior.noise.beta_0 = [0.01; 0.01; 0.01; 0.01];

prior.step_size = [0.003; 0.010; 0.030; 1.0; 1.5; ...
                   0.03; 0.05; 0.08; 0.08; 0.08; 0.06];

prior.mcmc.n_iterations = 5000;
prior.mcmc.burn_in = 500;
prior.mcmc.thin = 25;
prior.mcmc.n_chains = 4;
prior.mcmc.base_seed = 12345;




%prior.paths.base = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion';
prior.paths.base= '/Users/birotimi/Downloads/HVSR-Joint-Inversion';
prior.paths.cps_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_v3.30');
prior.paths.hv_dfa_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_v3.30_new');

end
