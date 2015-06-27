%%% Generate mu_kappa_psi plots %%%

delta = 0;
opts = struct();
opts.mod_psi = 0.6;
opts.max_iter = 100000;

mu_range = [-0.1,-1.3];
mu_points = 100;
mu_list = linspace(mu_range(1),mu_range(2),mu_points);

kappa_range = [0,0.1];
kappa_points = 100;
kappa_list = linspace(kappa_range(1),kappa_range(2),kappa_points);

mu_kappa_psi_data = struct();

%ALPHA = 0
alpha = [0,1];
lattice_dim = [1,2];
[M0, k_trans, k_max,diff_mat,err_mat] = mu_kappa_vs_psi(mu_list,kappa_list,...
            alpha,lattice_dim,delta,opts);
mu_kappa_psi_data('psi',M0,'mu_list',mu_list,'kappa_list',...
    kappa_list,'trans_list',k_trans,'max_list',k_max,...
    'diff_mat',diff_mat,'err_mat',err_mat);

save('./saves/mu_kappa_psi_0','mu_kappa_psi_data');
figure(1)
plot_mu_kappa_vs_psi(M0,mu_list,kappa_list,k_trans,k_max)

%ALPHA = 0.3
alpha = [3,10];
lattice_dim = [4,5];
[M3, k_trans, k_max,diff_mat,err_mat] = mu_kappa_vs_psi(mu_list,kappa_list,...
            alpha,lattice_dim,delta,opts);
mu_kappa_psi_data('psi',M3,'mu_list',mu_list,'kappa_list',...
    kappa_list,'trans_list',k_trans,'max_list',k_max,...
    'diff_mat',diff_mat,'err_mat',err_mat);
save('./saves/mu_kappa_psi_3','mu_kappa_psi_data');
figure(2)
plot_mu_kappa_vs_psi(M3,mu_list,kappa_list,k_trans,k_max)

%ALPHA = 0.6
alpha = [6,10];
lattice_dim = [4,5];
[M6, k_trans, k_max,diff_mat,err_mat] = mu_kappa_vs_psi(mu_list,kappa_list,...
            alpha,lattice_dim,delta,opts);
mu_kappa_psi_data('psi',M6,'mu_list',mu_list,'kappa_list',...
    kappa_list,'trans_list',k_trans,'max_list',k_max,...
    'diff_mat',diff_mat,'err_mat',err_mat);
save('./saves/mu_kappa_psi_6','mu_kappa_psi_data');
figure(3)
plot_mu_kappa_vs_psi(M5,mu_list,kappa_list,k_trans,k_max)

%ALPHA = 0.9
alpha = [9,10];
lattice_dim = [4,5];
[M9, k_trans, k_max,diff_mat,err_mat] = mu_kappa_vs_psi(mu_list,kappa_list,...
            alpha,lattice_dim,delta,opts);
mu_kappa_psi_data('psi',M9,'mu_list',mu_list,'kappa_list',...
    kappa_list,'trans_list',k_trans,'max_list',k_max,...
    'diff_mat',diff_mat,'err_mat',err_mat);
save('./saves/mu_kappa_psi_9','mu_kappa_psi_data');
figure(4)
plot_mu_kappa_vs_psi(M5,mu_list,kappa_list,k_trans,k_max)

M_diff = M5-M3;
figure(5)
plot_mu_kappa_vs_psi(M_diff,mu_list,kappa_list);
