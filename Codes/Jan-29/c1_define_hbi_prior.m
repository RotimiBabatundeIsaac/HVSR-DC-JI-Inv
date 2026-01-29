function prior = c1_define_hbi_prior()

prior.n_params = 19;
prior.param_names = {'h1','h2','h3','h4','h5','h6','h7','h8','h9',...
                     'Vs1','Vs2','Vs3','Vs4','Vs5','Vs6','Vs7','Vs8','Vs9','Vs10'};
prior.param_units = {'km','km','km','km','km','km','km','km','km',...
                     'km/s','km/s','km/s','km/s','km/s','km/s','km/s','km/s','km/s','km/s'};
prior.idx.h = 1:9;
prior.idx.vs = 10:19;
prior.enforce_monotonicity = true;
prior.monotonic_start_idx = 10;

h_min = [0.005; 0.02; 0.08; 0.30; 1.5; 2.0; 4.0; 6.0; 8.0];
h_max = [0.080; 0.30; 1.50; 4.00; 10.0; 15.0; 18.0; 22.0; 28.0];

vs_min = [0.12; 0.25; 0.50; 0.80; 2.00; 2.60; 3.00; 3.40; 3.80; 4.10];
vs_max = [0.90; 1.60; 2.20; 3.20; 3.80; 4.00; 4.30; 4.60; 4.90; 5.20];

prior.bounds.lower = [h_min; vs_min];
prior.bounds.upper = [h_max; vs_max];

prior.noise.n_noise = 4;
prior.noise.names = {'sigma_hvsr','sigma_ellip','sigma_cph','sigma_ugr'};
prior.noise.alpha_0 = [2.0; 2.0; 2.0; 2.0];
prior.noise.beta_0 = [0.01; 0.01; 0.01; 0.01];

prior.step_size = [0.003; 0.010; 0.030; 0.15; 0.40; 0.60; 0.80; 1.0; 1.2; ...
                   0.03; 0.05; 0.08; 0.10; 0.08; 0.07; 0.06; 0.06; 0.05; 0.05];

prior.mcmc.n_iterations = 5000;
prior.mcmc.burn_in = 500;
prior.mcmc.thin = 25;
prior.mcmc.n_chains = 4;
prior.mcmc.base_seed = 12345;

prior.paths.base = '/Users/birotimi/Downloads/HVSR-Joint-Inversion';
prior.paths.cps_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_v3.30');
prior.paths.hv_dfa_bin = fullfile(prior.paths.base, 'AGU_models', 'bin_v3.30_new');

end