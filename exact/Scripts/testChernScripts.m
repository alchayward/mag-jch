% Set up and run a Chern Number calculation
latticeDim = [4 4];
model = 'BH';
maxParticlesPerSite = 1;
nParticles = 2;
nAtoms=1;
alpha = [1 4];
onsiteStrength = [1,1];
hoppingStrength = 1;
ddStrength=0;
twist = [0 0];
A = [-1 0];
sysPath = './systems/';
atomLevels = 2;

degeneracy = 2;
nLevels = 3;
gridDims = [30 30];
fileName = 'test.mat';
filePath = './saves/test';
offSet = -15;

resume = 0;


system(sprintf('mkdir %s',filePath));
system(sprintf('cp ./ChernScript.m %sChernScripts.m',filePath));

sysParams=struct('model',model,'dim',latticeDim,'nParticles',nParticles,...
    'maxParticlesPerSite',maxParticlesPerSite,'nAtoms',nAtoms,...
    'atomLevels',atomLevels,'alpha',alpha, ...
    'onsiteStrength',onsiteStrength,...
    'ddStrength',ddStrength,'twist', twist,...
    'hoppingStrength',hoppingStrength,'A',[-1,0],'sysPath',sysPath);
    

chernParams = struct('degeneracy',degeneracy,'nLevels',nLevels,...
    'gridDims',gridDims,'fileName',fileName,'filePath',filePath,'offSet',...
offSet,'resume',resume);



h = HubbardLibrary;


grid = h.Chern(sysParams, chernParams);



