% Set up and run a Chern Number calculation


latticeDim = [5 4];
model = 'GJCH';
maxParticlesPerSite = 4;
nParticles = 4;
nAtoms=1;
alpha = [1,5];
onsiteStrength = 10*[1,sqrt(2),0,0]; 
hoppingStrength = 1;
ddStrength=0;
twist = [0 0];
A = [-1 0];
sysPath = './systems/';
atomLevels = 3;




degeneracy = 3;
nLevels = 3; %total levels = degen X nLevels


gridDims = [30 30];


fileName = 'B10D0v2.mat';
filePath = './saves/GJCH5x4p4Pfaff/';
offSet = -20;

resume = 0;

system(sprintf('mkdir %s',filePath));
system(sprintf('cp ./ChernScript.m %sChernScript.m',filePath));


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

%h.MakeSystem(sysParams,1);

grid = h.Chern(sysParams, chernParams);

fprintf('done with %s',fileName);

