function configs = Generate_Pfaffian_beta2_Configs(betas,lattice_configs,delta)

beta = 1;
kappa=1;

configs = cell([length(betas)*length(lattice_configs),1]);
ind = 1;
for ii = 1:length(lattice_configs)
   lc = lattice_configs{ii};
   
   %Create parameter configs. 
   for jj = 1:length(betas)
       dat = struct(); 
       dat.lattice_dims = lc.lattice_dims;
       dat.nParticles = lc.nParticles;
       dat.delta = delta;
       dat.beta_2 = betas(jj);
       dat.kappa = kappa;
       dat.beta = beta;
       dat.done=false;
       configs{ind} = dat;
       ind=ind+1;
   end
end