%%% Generate vortex plots %%%

delta = 0;
mu = -0.78; %tip of the first mott lobe
jc = meanfield_onsite_init(mu,delta);

kappa = 0.045;


opts = struct();
opts.mod_psi = 0.3;
opts.max_iter = 1000;


lattice_dims = [64,64];
alpha_q = prod(lattice_dims);

%Single Vortex

alpha_p = 1;
alpha = [alpha_p,alpha_q];
jcmf = meanfield_lattice_init(lattice_dims,alpha);
[ Psi1,diff,err ] = jch_lattice_meanfield(kappa,jcmf,jc,opts);

vortex_data = struct('Psi',Psi1,'alpha',alpha,...
    'mu',mu,'delta',delta,'diff',diff,'err',err); %#ok<NASGU>

save('./saves/vortex_1.mat','vortex_data');
figure(1)
psiPlot(Psi1,lattice_dims);

%Single Vortex

alpha_p = 64;
alpha = [alpha_p,alpha_q];
jcmf = meanfield_lattice_init(lattice_dims,alpha);
[ Psi64,diff,err ] = jch_lattice_meanfield(kappa,jcmf,jc,opts);

vortex_data = struct('Psi',Psi64,'alpha',alpha,...
    'mu',mu,'delta',delta,'diff',diff,'err',err);

save('./saves/vortex_64.mat','vortex_data');
figure(2)
psiPlot(Psi64,lattice_dims);