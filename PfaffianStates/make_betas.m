addpath('../exact/','../exact/systems')
nParticles = 4;
lattice = [4,4];
delta = -10;
betas = 0:.05:2; % about 40 points
beta_confs = Generate_Pfaffian_beta2_Configs(...
    betas,cellfun(@(ld)new_lattice_config(ld,nParticles),...
    {lattice},'UniformOutput',0),delta);
confs = Get_Config_Data(beta_confs);
filename = ['beta_data', sprintf('%d%d%d%d',nParticles,lattice(1),lattice(2),delta)];
save(filename,'confs');
