function [ jc] = meanfield_onsite_init( mu,delta)
%meanfield_iter_function returns the interation function for use in the
%meanfield convergence program.

    par=struct('onsiteStrength',struct('b',1,'D',delta,'mu',mu));

    jc=jch_onsite(par);
    jc.f_psi = jc.Make_fPsi_fitted();
end