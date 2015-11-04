addpath('../exact/','../exact/systems')
delta = -20;
beta_confs = Generate_Pfaffian_beta2_Configs(...
    0:.1:2,cellfun(@(ld)new_lattice_config(ld,4),...
    {[4,4],[4,5],[5,5],[5,6],[6,6]},'UniformOutput',0),delta);
confs = Get_Config_Data(beta_confs);
save('beta2_data','confs');
