% Set up and run a Chern Number calculation


lattice_dims = [6 6];
model = 'JCH';
maxParticlesPerSite = 4;
nParticles = 4;
nAtoms=1;
alpha = [8,36];
onsiteStrength = [1,0]; 
hoppingStrength = 1;
ddStrength=0;
twist = [0 0];
A = [-1 0];
sysPath = './systems/';
atomLevels = 2;




degeneracy = 2;
nLevels = 3; %total levels = degen X nLevels


gridDims = [30 30];


fileName = 'D0.mat';
filePath = './saves/JCH6x6p4Laughlin/';
offSet = -20;

resume = 0;

system(sprintf('mkdir %s',filePath));
system(sprintf('cp ./ChernScript.m %sChernScript.m',filePath));


sysParams=struct('model',model,'lattice_dim',lattice_dims,'nParticles',nParticles,...
    'maxParticlesPerSite',maxParticlesPerSite,'nAtoms',nAtoms,...
    'atomLevels',atomLevels,'alpha',alpha, ...
    'onsiteStrength',onsiteStrength,...
    'ddStrength',ddStrength,'twist', twist,...
    'hoppingStrength',hoppingStrength,'A',[-1,0],'sysPath',sysPath);
    

chernParams = struct('degeneracy',degeneracy,'nLevels',nLevels,...
    'gridDims',gridDims,'fileName',fileName,'filePath',filePath,'offSet',...
offSet,'resume',resume);



h = HubbardLibrary;

%h.MakeSystem(sysParams,1);

grid = h.Chern(sysParams, chernParams);

fprintf('done with %s',fileName);

