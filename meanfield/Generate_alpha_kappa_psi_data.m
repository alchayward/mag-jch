%%% Generate alpha_kappa plots %%%

delta = 0;
mu = -0.78; %tip of the first mott lobe
opts = struct();
opts.mod_psi = 0.6;
opts.max_iter = 100000;
alpha_range = [0,0.5]; % Only do half the range, because otherhalf is equal
alpha_points = 100;
alpha_p_list = linspace(alpha_range(1),alpha_range(2),alpha_points);

lattice_dims = [10,10];

kappa_range = [0,0.1];
kappa_points = 100;
kappa_list = linspace(kappa_range(1),kappa_range(2),kappa_points);



 [ psi,k_trans,k_max,diff_mat,err_mat] =...
     alpha_kappa_vs_psi(alpha_p_list,kappa_list,lattice_dims,delta,mu);
alpha_kappa_psi_data = struct('psi',psi,'alpha_list',alpha_list,...
    'kappa_list',kappa_list,'mu',mu,'delta',delta,...
    'trans_list',k_trans,'max_list',k_max,...
    'diff_mat',diff_mat,'err_mat',err_mat);

save('./saves/alpha_kappa_psi_tip','alpha_kappa_psi_data');
figure(1)
plot_mu_kappa_vs_psi(psi,alpha_p_list,kappa_list)