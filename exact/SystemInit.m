
p = struct(); %parameters structers
p.lattice_dim = [6 6];
p.model = 'GJCH';
p.maxParticlesPerSite = 4;
p.nParticles = 4;
p.nAtoms=1;
p.alpha = [4 36];
p.onsiteStrength = [1,sqrt(2),-5,-10];
p.hoppingStrength = 1;
p.ddStrength = 0;
p.twist = [0 0];
p.A = [-1 0];
p.sysPath = './systems/';
p.atomLevels = 3;

%nLevels = 6;
%filePath = './saves/testJCH.mat';
%strengths = [zeros(1,11);0:.1:1]';

%sysParams=struct('model',model,'dim',latticeDim,'nParticles',nParticles,...
%    'maxParticlesPerSite',maxParticlesPerSite,'nAtoms',nAtoms,...
%    'alpha',alpha, 'onsiteStrength',onsiteStrength,'twist', twist,...
%    'hoppingStrength',hoppingStrength,'ddStrength',ddStrength,...
%    'A',[-1,0],'sysPath',sysPath);
    
h=HubbardLibrary;
s=h.Lattice(p);
s = s.Initilize();
