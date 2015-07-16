% Set up and run a Chern Number calculation

latticeDim = [4 4];
model = 'BH';
maxParticlesPerSite = 4;
nParticles = 4;
nAtoms=2;
alpha = [1 4];
hoppingStrength = 1;
twist = [0 0];
A = [-1 0];
sysPath = './systems/';


degeneracy = 3;
nLevels = 3;
gridDims = [30 30];
filePath = './saves/BH4x4p4Pfaff/list/';
offSet = -15;
stConf={};
fnConf={};
for ii = 1:11
for jj = 1:11
   
   stConf{ii,jj} = [0,(ii-1)*.1,(jj-1)*.1];
   fnConf{ii,jj} = ['bhu', num2str(ii),'v',num2str(jj),'.mat'];
   

end
end

h = HubbardLibrary;
grid={};
for ii = 1:11
for jj = 1:11
 onsiteStrength = stConf{ii,jj};
 fileName = fnConf{ii,jj};

sysParams=struct('model',model,'dim',latticeDim,'nParticles',nParticles,...
    'maxParticlesPerSite',maxParticlesPerSite,'nAtoms',nAtoms,...
    'alpha',alpha, 'onsiteStrength',onsiteStrength,'twist', twist,...
    'hoppingStrength',hoppingStrength,'A',[-1,0],'sysPath',sysPath);
    

chernParams = struct('degeneracy',degeneracy,'nLevels',nLevels,...
    'gridDims',gridDims,'fileName',fileName,'filePath',filePath,'offSet',...
offSet);
disp(mat2str([ii,jj]))
grid{ii,jj} = h.Chern(sysParams, chernParams);
end
end
