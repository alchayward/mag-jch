detunings = -10:1:10;
lattice = [4,4]; 
nParticles = 2;
addpath('../exact/','../exact/systems')
de_confs = Generate_Laughlin_Detuning_Configs(...
    detunings,cellfun(@(ld)new_lattice_config(ld,nParticles),...
    {lattice},'UniformOutput',0));
confs = Get_Config_Data(de_confs);
filename = ['delta_data', arrayfun(@(x)num2str(x),[lattice(1),lattice(2),nParticles])];
save(filename,'confs');
