addpath('../exact/','../exact/systems')
nParticles = 4;
lattice = [4,4];
delta = -20;
beta_confs = Generate_Pfaffian_beta2_Configs(...
    0:.1:2,cellfun(@(ld)new_lattice_config(ld,4),...
    {[4,4],[4,5],[5,5],[5,6],[6,6]},'UniformOutput',0),delta);
confs = Get_Config_Data(beta_confs);
filename = ['beta_data', arrayfun(@(x)num2str(x),[lattice(1),lattice(2),nParticles,delta])];
save(filename,'confs');
