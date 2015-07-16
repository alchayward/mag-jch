
latticeDim = [6 6];
model = 'BH';
maxParticlesPerSite = 2;
nParticles = 4;
nAtoms=1;
alpha = [4 36];
onsiteStrength = [0,0];
hoppingStrength = 1;
ddStrength = 0;
twist = [0 0];
A = [-1 0];
sysPath = './systems/';

%nLevels = 6;
%filePath = './saves/testJCH.mat';
%strengths = [zeros(1,11);0:.1:1]';

sysParams=struct('model',model,'dim',latticeDim,'nParticles',nParticles,...
    'maxParticlesPerSite',maxParticlesPerSite,'nAtoms',nAtoms,...
    'alpha',alpha, 'onsiteStrength',onsiteStrength,'twist', twist,...
    'hoppingStrength',hoppingStrength,'ddStrength',ddStrength,...
    'A',[-1,0],'sysPath',sysPath);
    
h=HubbardLibrary;
sb2=h.Lattice(sysParams);
