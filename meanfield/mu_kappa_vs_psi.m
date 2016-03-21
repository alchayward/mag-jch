function [ psi,kappa_t,kappa_max,diff_mat,err_mat ] =...
    mu_kappa_vs_psi(mu_list,kappa_list,...
    alpha,lattice_dims,delta,varargin)
%mu_kappa_vs_psi 
%mu_list: list of values of the chemical to solve for.
%kappa_list: list of values of the hopping strength to solve for.
%alpha: [p,q], a list of intergers, q>0, s.t. alpha = p/q.
%lattice_dims: lattice dimensions [lx,ly].
%delta: detuning parameter.
%opts: options that get sent to the solver. see ConvergePsi
opts = struct();
    if nargin > 5
        opts = varargin{1};
    end


psi = zeros(length(mu_list),length(kappa_list));
kappa_t = 0*mu_list;
kappa_max = 0*mu_list;
diff_mat = psi*0;
err_mat = psi*0;


jcmf = meanfield_lattice_init(lattice_dims,alpha);
%matlabpool(2) %parallel is not working on my comp right now.
for ii = 1:length(mu_list)
    mu = mu_list(ii);
    jc = meanfield_onsite_init(mu,delta);
    fprintf('mu = %d\n',mu)
    
    [psi(ii,:),kappa_t(ii),kappa_max(ii),diff_mat(ii,:),err_mat(ii,:)] =...
        psi_vs_kappa(kappa_list,jcmf,jc,opts);
    
    
end
%matlabpool close
end

