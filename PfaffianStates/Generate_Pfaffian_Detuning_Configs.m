function configs = Generate_Pfaffian_Detuning_Configs(detunings,lattice_configs)
beta_2=sqrt(2);
beta = 1;
kappa=1;

configs = cell([length(detunings)*length(lattice_configs),1]);
ind = 1;
for ii = 1:length(lattice_configs)
   lc = lattice_configs{ii};
   
   %Create parameter configs. 
   for jj = 1:length(detunings)
       dat = struct(); 
       dat.lattice_dims = lc.lattice_dims;
       dat.nParticles = lc.nParticles;
       dat.delta = detunings(jj);
       dat.beta_2 = beta_2;
       dat.kappa = kappa;
       dat.beta = beta;
       dat.done=false;
       configs{ind} = dat;
       ind=ind+1;
   end
end

