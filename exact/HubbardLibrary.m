function h = HubbardLibrary()
%HubbardLibrary: Library of functions for Hubbard Model Calculation
%   Functions: Chern(systemParameters, chernParameters)
%              
h.Chern = @Chern;
h.Lattice = @Lattice;
h.ChernCalc = @ChernCalc;
h.GetHamiltonianAndTwistMatricies = @GetHamiltonianAndTwistMatricies;
h.MakeSystem = @MakeSystem;
h.SortEigs = @SortEigs;
h.VaryOnsite = @VaryOnsite;
h.GetSys = @GetSys;
h.GramSchmit = @GramSchmit;
h.ArnoldiVectors = @ArnoldiVectors;
h.Groundstate = @Groundstate;
end


function grid = Chern( systemParameters,chernParameters)
%CHERN Summary of this function goes here
%   Detailed explanation goes here

%check that directory exists, if not create it.  
if 7 ~= exist(chernParameters.filePath,'file')
    mkdir(chernParameters.filePath);  
end

opFile = fullfile(chernParameters.filePath,'./chernOperators.mat');
saveFile = fullfile(chernParameters.filePath,...
            chernParameters.fileName);
tempFile = fullfile(chernParameters.filePath,...
            'tempChern.mat');
outFile = fullfile(chernParameters.filePath,...
            'chern.out');

[h, tx, ty] = ...
    GetHamiltonianAndTwistMatricies(systemParameters); %#ok<*NASGU,*ASGLU>
save(opFile,'h','tx','ty');
clear('h','tx','ty');
gds = chernParameters.gridDims;
gridPositions = zeros(prod(gds),2);

%This is dumb code. It will work.
dth = 2*pi./gds;
[x,y] = meshgrid(dth(1):dth(1):2*pi,dth(2):dth(2):2*pi);

gridPositions(:,1) = reshape(x,prod(gds),1);
gridPositions(:,2) = reshape(y,prod(gds),1);

refStatePos =[0,0;pi,pi];
   
grid = ChernCalc(opFile, gridPositions, refStatePos, ...
    chernParameters.degeneracy,chernParameters.nLevels,...
    tempFile,outFile,chernParameters.offSet,chernParameters.resume);
save(saveFile,'grid','systemParameters','chernParameters');
end

function obj = Lattice(params)
% params should include:
%struct('model','dim','pn','alpha',
%    'onsitestrength','hoppingstrength','twist','A',[-1,0]);

switch params.model
	case 'BH'
		obj = BHsystem(params);
    case 'BHLR'
        obj = BHLRsystem(params);
	case 'GJCH'
		obj = GJCHsystem(params);
    case 'JCH'
		obj = JCHsystem(params);
	case 'HC'
		obj = HCsystem(params);
    case 'TCH'
		obj = TCHsystem(params);
end

end
                     
function grid = ChernCalc(operatorFile, gridPositions,...
                referenceStatePositions, groundStateDegeneracy,...
                nLevels,tempFile,outFile,offSet,varargin)
         
%%%                     CHERNGRID
%%%
%%%
%%%
%%% Computes the auxilory fuctions defined in [reference] for the
%%% computation of chern numbtrs associated with a hamiltonian. 
%%%
%%%
%%%
%%%
%%% function grid = CHERNGRID(operatorFile, gridPositions,...
%%%               referenceStatePositions, groundStateDegeneracy,...
%%%                nLevels,tempFile,varargin)
%%%
%%%                         INPUTS
%%%
%%%
%%%  operatorFile(string): Filepath for the hamiltonian and twist
%%%  matricies.
%%%  gridPostions(Mx2 double): Array with the m-th entry being x and y
%%%  twist positions.  
%%%  
%%%     energyLowerBound
%%%    
%%%  referenceStatePostions(2x2 double): a 2*2 array with the two reference
%%%  state twist
%%%  angles. 
%%%  
%%%  groundStateDegeneracy(int): degeneracy of the ground state manifold.      
%%%
%%%  nLevels(int): How many levels to compute. The total number of
%%%  eignestates found will then be groundStateDegeneracy*nLevels.
%%%
%%%  tempFile(string): filepath to temporary store the results
%%%  of the computation.
%%%
%%%  resume(bool)(OPTIONAL): 
%%%  
%%%  iterationPartitionMultiplier(int)(OPTIONAL): 
%%%
%%%
%%%
%%%                         OUTPUTS

    outFileID = 2;
    eigTriesMax = 100;
    if outFile ~= 2
        outFileID = fopen(outFile,'w');
    end
    
    global h tx ty
    fprintf(outFileID,'Begin Chern Number Calc\n');
    
    % only want 2 optional inputs at most
    numvarargs = length(varargin);
    if numvarargs > 2
        error('myfuns:somefun2Alt:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    %%% Set Defaults for optional inputs
    resume = 0;
    iterationPartitionMultiplier = 30;
    numberOfWorkers=1;
    if numvarargs > 0
        resume = varargin{1};
    end
    if numvarargs > 1
        iterationPartitionMultiplier = varargin{2};
    end
    q = groundStateDegeneracy;

    %%% eigs paramters:
    %number of eigenstates to compute.        
    nEigenstates = groundStateDegeneracy*nLevels;
    fprintf(outFileID,'Load Matricies\n');       
    %%%Get Hamiltonian and Twist  Matricies
    ops = load(operatorFile);
    h = ops.h;
    tx = ops.tx;
    ty = ops.ty;
    
    clear('ops');
    fprintf(outFileID,'Matricies Loaded\n');    

    %%%Compute Reference States
    hxy = MakeTwistedHamiltonian(referenceStatePositions(1,:));
    %First get an estimate for the groundstates with eigs to give to
    %lobpcg)
    fprintf(outFileID,'Compute First Reference Groundstate\n');   
    
    %[V D] = eigs(hxy,nEigenstates,'SR');
    %[V ~] = SortEigs(V,diag(D));
     %%%trying my implimentation ofarnodi iteration. eigs takes too long.
     aOpts=struct('nLevels',nEigenstates,'nIterations',100,'offSet',offSet);
     [ref1, D, ~] = ArnoldiVectors(h,aOpts);
     fprintf(outFileID,'Compute First Reference Groundstate using lobpcg\n');   
     lhist=0;
     rnhist=0;
     [ref1,D , failureFlag]=lobpcg(ref1,hxy);   
     eigTries = 0;
     while failureFlag == 1 && eigTries < eigTriesMax   
        [ref1, D, ~] = ArnoldiVectors(hxy,aOpts,ref1);
        [ref1, D, failureFlag]=lobpcg(ref1,hxy);
        eigTries = eigTries +1;
    end
     %%% Need to do something about checking for convergence. If it fails, is
     %there a fallback? 
     hxy = MakeTwistedHamiltonian(referenceStatePositions(2,:));
     fprintf(outFileID,'Compute Second Reference Groundstate using lobpcg\n'); 
     [ref2, D, failureFlag]=lobpcg(ref1,hxy);    
    eigTries = 0;
    while failureFlag == 1 && eigTries < eigTriesMax  
        [ref2, D, ~] = ArnoldiVectors(hxy,aOpts,ref2);
        [ref2, D, failureFlag]=lobpcg(ref2 ,hxy);  
        eigTries = eigTries +1;
    end
    
     
   
    
%%%Initiate Grid 
     gridDim = size(gridPositions,1);
     E = zeros(gridDim,nEigenstates);
     th = zeros(gridDim,nLevels);
     A1 = zeros(gridDim,nLevels);
     A2 = zeros(gridDim,nLevels);
     success = zeros(gridDim,nLevels);
    %%%Compute grid
    
    %Set ground state starting vector:
    %eigsopts.v0 = ref1(:,1);
    
    %%% Set up loop so it can intermitantly save the computation, and so it
    %%% can be resumed midway if the function fails to complete for some
    %%% reason. The number of iterations per partition has to at least be
    %%% greater than the number of matpool workers. 
    
    %%% Work out where to start from
    if resume == 1
        T = load(tempFile);
        E = T.E;
        th = T.th;
        A1 = T.A1;
        A2 = T.A2;
        ref1 = T.ref1;
        ref2=T.ref2;
        startFrom = T.ii+1;
        clear('T');  
    else
        startFrom = 0; 
    end
    
    iterationsPerPartition = iterationPartitionMultiplier*numberOfWorkers;
    numberOfPartitions = ceil(...
        (gridDim-startFrom+1)/iterationsPerPartition);
    
    MTHhandle = @MakeTwistedHamiltonian;
    tStart=tic; %%%BEGIN THE COUNTDOWN!
    
    mPoolSize = matlabpool('size');
    if mPoolSize ~= 0
        matlabpool('close','force','local');
    end        
    matlabpool('open');
    mPoolSize = matlabpool('size');
    
    fprintf(outFileID,'matlab pool size: %d\n',mPoolSize);
    fprintf(outFileID,'Start Loop\n'); 
    for ii = startFrom:(numberOfPartitions-1)
        loopstart = 1 + ii*iterationsPerPartition;
        loopfinish = min([(ii+1)*iterationsPerPartition,gridDim]);
        
        %Loop over twist angles.
        parfor jj = loopstart:loopfinish
            %Compute h(x,y). Multiply h by x and y twist matricies.
            hxy = feval(MTHhandle,gridPositions(jj,:));
            
            [V, E(jj,:), failureFlag]=lobpcg(ref1,hxy); %Compute the Eigenstates
            eigTries = 0;
            while failureFlag == 1 && eigTries < eigTriesMax
                [V, D, ~] = ArnoldiVectors(hxy,aOpts,V);
                [V, D, failureFlag]=lobpcg(V,hxy);
                eigTries = eigTries +1;
            end
            %%%Compute the auzilory functions for each eigenstate manifold.
            
            s1 = V'*ref1;  %overlaps for each
            s2 = V'*ref2;  %reference state.
            %Loop over manifolds
            thp = zeros(nLevels,1);
            A1p = zeros(nLevels,1);
            A2p = zeros(nLevels,1);
            for kk = 1:nLevels;
                s1n=s1((q*(kk-1) + 1):(q*(kk)),(q*(kk-1) + 1):(q*(kk)));
                s2n=s2((q*(kk-1) + 1):(q*(kk)),(q*(kk-1) + 1):(q*(kk)));
                thp(kk)=det(s1n'*s2n);
                thp(kk)=thp(kk)./abs(thp(kk));
                A1p(kk)=det(s1n'*s1n);
                A2p(kk)=det(s2n'*s2n);
            end
                th(jj,:)=thp;
                A1(jj,:)=A1p;
                A2(jj,:)=A2p;
        end
        tElapsed=toc(tStart);
        save(tempFile,'gridPositions','E','th','A1','A2','ref1','ref2',...
                    'loopfinish','ii');
        fprintf(outFileID,strcat('estimated time left: ',...
            num2str((gridDim - loopfinish)*tElapsed/(60^2*loopfinish))...
            ,' hours\n'));
    end
    fprintf(outFileID,'DONE!!!!!!!!!!!!!!!!!!!!!!\n');
    grid = struct('x',gridPositions(:,1),'y',gridPositions(:,2),...
            'th',th,'A1',A1,'A2',A2,'E',E);
    matlabpool close
    fclose(outFileID);
    
    function hxy =  MakeTwistedHamiltonian(angle)    
        %So what this does is pick out the twisted hopping elements and
        %multiplies them by the twist phase, for which we have to do
        %complex conjugates seperately. 
        tax = exp(1i*angle(1));
        tay = exp(1i*angle(2));
        hxy = h;
        hxy(tx) = hxy(tx).*tax;
        hxy(tx') = hxy(tx').*conj(tax);
        hxy(ty) = hxy(ty).*tay;
        hxy(ty') = hxy(ty').*conj(tay);
    end
end

function [h, tx, ty] = ...
    GetHamiltonianAndTwistMatricies(systemParameters,varargin)



%%%
%%%
%%%
%%%
%%%
%%%                     GetHamiltonianAndTwistMatricies
%%%
%%%
%%%
%%% Generate the hamiltonain and twist matricies needed for the chern
%%% number calculation.
%%%
%%%
%%%
%%% function [h txp txm typ tym] = ...
%%%   GetHamiltonianAndTwistMatricies(systemParameters,newSystem)
%%%
%%%
%%%
%%%
%%%                         INPUTS
%%%  systemParameters(struct): a struct containing the following
%%%  variables:
%%%
%%%     model(string):
%%%     dim(2 x int ): x and y dimesions of the lattice
%%%     nParticles(int): number of excitations in the lattice
%%%     alpha(2 x int): alpha is the flux paramters = p/q
%%%     onsiteStrength(vector double): array of onsite energies
%%%     hoppingStrength(double): inter-site hubbard hopping energy.
%%%     twist(DEPRICATED) (2x double): twist angle of the hopping matrix
%%%     gauge(2*double): Gauge of the vector potential. ie
%%%     A=alpha*(gauge(1) y , gauge(2) x ). Should default to [-1,0];
%%%     path(string): location of storage directory.
%%%
%%%  newSystem(bool)(OPTIONAL default = 0): Generate a system from scratch,
%%%                  even if it already exists.
%%%
%%%                         OUTPUTS
%%%
%%%  h(MxM c-double): hamiltonain matrix.
%%%  txp(MxM double): 4 twist matricies.
%%%
%%%
%%%



% only want 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional inputs');
end


%%% Set Defaults for optional inputs


%sysPath = systemParameters.sysPath;

file = MakeSystem(systemParameters,0);

sys = load(file,'sys');  % load system from file.
sys = sys.sys; %A HACK

%%% Get the matricies using the system class functions.
h = sys.MakeHamiltonian(systemParameters.hoppingStrength,...
    systemParameters.onsiteStrength);
[tx, ty] = sys.MakeTwistMatricies;
clear('sys');
end

function file = MakeSystem(systemParameters,varargin)
%varargin is a bool, if 1, make a new matricies even if file
%already exisits.
newMats = 0;
if length(varargin) == 1
    newMats = varargin{1} ;
end
sys = Lattice(systemParameters);
file = sys.GenerateFilePath;
if (exist(file,'file') ~= 2) || newMats
    disp('Making New System');
    sys = sys.Initilize;
    save(file,'sys');
end
end

function [V, D] = SortEigs(V, d)
[D, ix] = sort(real(d));
V=V(:,ix);

V = GramSchmit(V); %Orthonormilize
end


function V = GramSchmit(V)
V(:,1) = V(:,1)/sqrt(V(:,1)'*V(:,1));
for ii = 2:size(V,2);
    for jj = 1:(ii-1)      
        proj = V(:,jj)'*V(:,ii);
        V(:,ii) = V(:,ii) - proj*V(:,jj);
        V(:,ii) = V(:,ii)/sqrt(V(:,ii)'*V(:,ii));
    end
end
end

function E = VaryOnsite(systemParameters, strengths, nLevels, saveFilePath,outFilePath, varargin )

%		     VaryOnsite
%
% function E = VaryOnsite(systemParameters, strengths,saveFileName, varargin )
%
%

eigTriesMax = 100;

% only want 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

    offSet=0;
if numvarargs >= 1
    offSet = varargin{1};
end

nOperators = 1;
whichOperator = 'h';



outFileID = 2;

if outFilePath ~= 2
    outFileID = fopen(outFilePath,'w');
end



if numvarargs >= 2
    operators = varargin{2};
    nOperators = length(operators);
    whichOperator = 'm';
end



%nOperators = length(operators);


nStrengths=size(strengths,1);

%%% Get the hamiltonian using the system class functions.

sys=GetSys(systemParameters);

E = zeros(nStrengths,nLevels,nOperators);

fprintf(outFileID,'Start Loop\n'); 
for ii = 1:nStrengths
    
     fprintf(outFileID,'loop iteration number: %d\n',ii); 
    
    h = sys.MakeHamiltonian(systemParameters.hoppingStrength,...
        strengths(ii,:));
    
    
    if 1 == ii
        fprintf(outFileID,'First Eigenstate\n'); 
        %[V D] = eigs(h,nLevels,'SR');
        %[V D] = SortEigs(V, D);
        
        aOpts=struct('nLevels',nLevels,'nIterations',50,'offSet',offSet);
        
        [V, D, ~] = ArnoldiVectors(h,aOpts);
        [V, D, failureflag]=lobpcg(V,h);
        
        fprintf(outFileID,'First found\n'); 
    end
    
    failureFlag = 1;
    eigTries = 0;
    while failureFlag == 1 || eigTries < eigTriesMax  
    
        [V, D, failureFlag]=lobpcg(V,h);
        eigTries = eigTries +1;
    end
    
    if failureFlag == 1
        %fprint(outFileID,('Convergence failed for str:' + num2str()
    end
    if whichOperator == 'h'
        E(ii,:) = D;
    else
        for jj = 1:nOperators
            E(ii,:,jj) = V'*operators{jj}*V;
        end
    end
    
end

 fprintf(outFileID,'Calc Done!\n');

save(saveFilePath,'E','strengths','systemParameters');
fclose(outFileID);

end

function [V, D, flag] = Groundstate(H,opts,varargin)
%Find the lowest eigenstates of H.
%Several methods are availiable by setting the method option:
%   'eigs' : use the matlab function eigs(A,k,'SR') (default)
%    'lobpcg' : gradient method
%    'arnoldi' : use my defintion of the gram-schmit arnoldi method
%    'mix' : uses a mixture of lobpcg and aroldi, which tries to avoid the
%           lobpcg method from getting stuck
% other options 
%     'k(8)' : number of levels
%     'iterations(50)'
%     'accuracy'
%     'V0'(random) : initial trial groundstate. must be length(A) by k
%     unless eigs is being used, in which case, it is length(A) by 1.
%     



if ~isfield(opts,'k')
    opts.k = 6;
end
    nLevels = opts.k;

if ~isfield(opts,'tol')
    opts.tol = 10^-12;
end
tol = opts.tol;

if ~isfield(opts,'nIterations')
    opts.nIterations = 50;
end
nIterations = opts.nIterations;

if isfield(opts,'offSet')
    H = H + opts.offSet*speye(size(A,1));
    offset = opts.offset;
else
    opts.offset = 0;
    offset = opts.offset;
end

if ~isfield(opts,'method')
    opts.method = 'eigs';
    method = 'eigs';
end
method = opts.method;

if ~isfield(opts,'V0')
    if strcmp(method,'eigs')
        V0n = 1;
    else
        V0n = nLevels;
    end
    opts.V0 = rand(length(H),V0n);
end
V0 = opts.V0;

if strcmp(method,'eigs')
    aOpts = struct('V0',V0,'tol',tol);
    try
        [V, D, flag] = eigs(H,nLevels,'SR',aOpts);
        D=diag(D);
    catch ME
       disp('eigs failed to converge, trying another');
       new_opts = opts;
       new_opts.method = 'lobpcg';
       [V, D, flag] = Groundstate(H,new_opts,varargin);
    end
       
elseif strcmp(method,'lobpcg')
    flag = 1;
    eigTries = 0;
    failureFlag = 1;
    while failureFlag == 1 && eigTries < nIterations  
        [V0, D, flag]=lobpcg(V0,H);
        failureFlag=flag;
        eigTries = eigTries +1;
    end
V = V0;
elseif strcmp(method,'arnoldi') 
        aOpts=struct('nLevels',k,'nIterations',50,'offSet',offSet);
        [V, D, flag] = ArnoldiVectors(h,aOpts);
elseif strcmp(method,'mix')
end
V = GramSchmit(V);
D=D-offset;
%%% sort out some flag stuff (failure to converge, ect.)

%%% Sort Eigenvectors into order (sometimes they are not in order)
[D,order] = sort(D,'ascend');
V = V(:,order);
end

function sys = GetSys(systemParameters)
file = MakeSystem(systemParameters);
sys = load(file,'sys');  % load system from file.
sys = sys.sys; %A HACK
end


function [V, D, diff] = ArnoldiVectors(A,opts,varargin)
%%% Returns the largest magnitude 
    if isfield(opts,'nLevels')
        nLevels = opts.nLevels;
    else
        nLevels = 8;
    end
    % 
    % if isfield(opts,'V0')
    % 
    if isfield(opts,'tol')
        tol = opts.tol;
    else
        tol = 10^-12;
    end

    if isfield(opts,'nIterations')
        nIterations = opts.nIterations;
    else
        nIterations = 50;
    end

    if isfield(opts,'offSet')
        A = A + opts.offSet*speye(size(A,1));
    end

    if ~isempty(varargin)
        V = varargin{1};
    else
        V = rand(size(A,1),nLevels);
    end

    V = GramSchmit(V);
    D=zeros(nLevels,1);
    for ii = 1:nIterations
        Vn = A*V;
        d = diag(V'*Vn);
        diff = max(abs(d-D));
        if diff < tol
            break;
        end
        [V, D] = SortEigs(Vn,d);
    end

    if isfield(opts,'offset')
        D = D - opts.offSet;
    end
end
