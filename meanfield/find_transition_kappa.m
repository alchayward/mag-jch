function [ kappa ] = find_transition_kappa(varargin)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    f_psi = varargin{1};
    jcmf = varargin{2};
    amat = jcmf.amat;
elseif nargin == 1
    %Take a structure with all the parameters (mu,alpha_p,alpha_q,delta)
    params = varargin{1};
    par_onsite=struct('onsiteStrength',struct('b',1,'D',params.delta,...
        'mu',params.mu));
    jc=jch_onsite(par_onsite);
    f_psi = jc.Make_fPsi();
    
    shift = 0;
    par_lattice=struct('dim',params.lattice_dims,'alpha',...
        params.alpha,'yshift',shift);
    jcmf=jchmf(par_lattice);
    amat = jcmf.amat;
elseif nargin == 4
    delta = varargin{1};
    mu = varargin{2};
    alpha = varargin{3};
    lattice_dims = varargin{4};
    par_onsite=struct('onsiteStrength',struct('b',1,'D',delta,...
        'mu',mu));
    jc=jch_onsite(par_onsite);
    f_psi = jc.Make_fPsi();
    
    shift = 0;
    par_lattice=struct('dim',lattice_dims,'alpha',...
        alpha,'yshift',shift);
    jcmf=jchmf(par_lattice);
    amat = jcmf.amat;
    
else
    msgID = 'find_transition_kappa:BadParamters';
    msg = 'Wrong number of parameters';
    baseException = MException(msgID,msg);
    throw(baseException)
end

%Compute the f_psi gradient.
accuracy = 10^-8;


a = f_psi(accuracy)/accuracy;

%This seems to fail sometimes. Can reduce min error. or catch exceptions
opts = struct('v0',jcmf.psi0);
try
    e = max(abs(eigs(amat,1,'lm',opts)));
catch ME
    opts.tol = 10^-4;
    try 
        e = max(abs(eigs(amat,1,'lm',opts)));
    catch ME2
        e = 4;
    end
end

kappa = 1/(e*a);

end

