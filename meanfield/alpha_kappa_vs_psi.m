function [ psi,kappa_t,kappa_max,diff_mat,err_mat] =...
    alpha_kappa_vs_psi(alpha_p_list,kappa_list,lattice_dims,delta,mu,opts)
%alpha_kappa_vs_psi 
%alpha_p_list: list of integers for alpha_p. alpha_q is the number of
%lattice sites. (so alpha_p = number of vorticies).
%kappa_list: list of values of the hopping strength to solve for.
%lattice_dims: lattice dimensions [lx,ly].
%delta: detuning parameter.
%mu: value of the chemical potential to solve for.
psi = zeros(length(alpha_p_list),length(kappa_list));



jc = meanfield_onsite_init(mu,delta);
alpha_q = prod(lattice_dims);

kappa_t = alpha_p_list*0;
kappa_max = alpha_p_list*0;
diff_mat = psi*0;
err_mat = psi*0;

parpool(4)
for ii = 1:length(alpha_p_list)
    alpha = [alpha_p_list(ii),alpha_q];
    jcmf = meanfield_lattice_init(lattice_dims,alpha);

    fprintf('alpha = %d\n',alpha(1)/alpha(2))
    [psi(ii,:),kappa_t(ii),kappa_max(ii),diff_mat(ii,:),err_mat(ii,:)]=... 
    psi_vs_kappa(kappa_list,jcmf,jc,opts);
    
end
parpool close
end
