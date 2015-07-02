function [ psi,diff,err ] = jch_lattice_meanfield(kappa,jcmf,jc,varargin)
%jch_lattice_meanfield: solve the mean field of a jch lattice.
%WARNING! There is nothing in here that checks for the divirgent criteria
    %( z*kappa > mu). so yeah, be careful.
    opts = struct();
    if nargin > 3
        opts = varargin{1};
    end
    
    f_psi_kappa = @(x)jc.f_psi(kappa*jcmf.amat*x);
    opts.lattice_dims = jcmf.lattice_dims;
    
    if ~isfield(opts,'psi0')
        opts.psi0 = rand(prod(jcmf.lattice_dims),1)*10^-6;
    end
    psi0 = opts.psi0;
    
    [psi,diff,err] = ConvergePsi(psi0,f_psi_kappa,opts);
    %plotPsi(psi,jcmf.latticeDim);
end