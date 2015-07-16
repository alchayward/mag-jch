% Set up and run a Vary strength calculation

latticeDim = [4 4];
model = 'BH';
maxParticlesPerSite = 4;
nParticles = 4;
nAtoms=2;
alpha = [1 4];
onsiteStrength = [0,0];
hoppingStrength = 1;
twist = [pi pi];
A = [-1 0];
sysPath = './systems/';


sysParams=struct('model',model,'dim',latticeDim,'nParticles',nParticles,...
    'maxParticlesPerSite',maxParticlesPerSite,'nAtoms',nAtoms,...
    'alpha',alpha, 'onsiteStrength',onsiteStrength,'twist', twist,...
    'hoppingStrength',hoppingStrength,'A',[-1,0],'sysPath',sysPath);
 

nLevels = 6;
filePath = './saves/BH4x4p4Pfaff/Vary1.mat';
outPath =  './saves/BH4x4p4Pfaff/vary1.o';

vars = 0:.01:1;
o = ones(1,length(vars));
strengths = [0*o;0*o;vars]';
offSet = -20;
   

%varyParams = struct('nLevels',nLevels,...
%    'fileName',fileName,'filePath',filePath);



h = HubbardLibrary;

E = h.VaryOnsite(sysParams, strengths, nLevels, filePath,outPath,offSet);
