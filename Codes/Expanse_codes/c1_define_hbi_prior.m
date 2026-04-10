function prior = c1_define_hbi_prior_10layer()

prior.n_params = 19;

prior.param_names = {'h1','h2','h3','h4','h5','h6','h7','h8','h9', ...
                     'Vs1','Vs2','Vs3','Vs4','Vs5','Vs6','Vs7','Vs8','Vs9','Vs10'};

prior.param_units = {'km','km','km','km','km','km','km','km','km', ...
                     'km/s','km/s','km/s','km/s','km/s','km/s','km/s','km/s','km/s','km/s'};

prior.idx.h  = 1:9;
prior.idx.vs = 10:19;

prior.enforce_monotonicity = true;
prior.monotonic_start_idx  = 10;

h_min = [0.005; 0.02; 0.08; 0.30; 1.50; 2.00; 4.00; 6.00; 8.00];
h_max = [0.080; 0.30; 1.50; 4.00; 10.0; 15.0; 18.0; 22.0; 28.0];

vs_min = [0.12; 0.25; 0.50; 0.80; 2.00; 2.60; 3.00; 3.40; 3.80; 4.10];
vs_max = [0.90; 1.60; 2.20; 3.20; 3.80; 4.00; 4.30; 4.60; 4.90; 5.20];

prior.bounds.lower = [h_min;  vs_min];
prior.bounds.upper = [h_max;  vs_max];

prior.depth_min_km = sum(h_min);
prior.depth_max_km = sum(h_max);

prior.noise.n_noise = 4;
prior.noise.names   = {'sigma_hvsr','sigma_ellip','sigma_cph','sigma_ugr'};
prior.noise.alpha_0 = [2.0; 2.0; 2.0; 2.0];
prior.noise.beta_0  = [0.01; 0.01; 0.01; 0.01];

prior.noise.use_correlated = true;
prior.noise.corrlen_names = {'L_hvsr','L_ellip','L_cph','L_ugr'};
prior.noise.corrlen_lb    = [0.01; 0.01; 0.01; 0.01];
prior.noise.corrlen_ub    = [1.50; 1.50; 1.50; 1.50];
prior.noise.corrlen_step  = [0.05; 0.05; 0.05; 0.05];
prior.noise.corrlen_init  = [0.10; 0.10; 0.10; 0.10];

prior.step_size = [ ...
    0.003; 0.010; 0.030; 0.15; 0.40; 0.50; 0.60; 0.80; 0.80; ...
    0.02;  0.03;  0.05;  0.06; 0.06; 0.05; 0.05; 0.05; 0.05; 0.05  ...
];

prior.mcmc.n_iterations = 5000;
prior.mcmc.burn_in      = 500;
prior.mcmc.thin         = 25;
prior.mcmc.n_chains     = 4;
prior.mcmc.base_seed    = 12345;

prior.paths.base      = '/expanse/lustre/projects/syu127/brotimi/HVSR-Joint-Inversion';
prior.paths.cps_bin   = '/expanse/lustre/projects/syu127/brotimi/Inversion/bin';
prior.paths.hv_dfa_bin= '/expanse/lustre/projects/syu127/brotimi/Inversion/bin';

end
