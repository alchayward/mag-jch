function [ jcmf] = meanfield_lattice_init(lattice_dims,alpha)
%meanfield_iter_function returns the interation function for use in the
%meanfield convergence program.
    shift = 0;
    par=struct('lattice_dims',lattice_dims,'alpha',alpha,'yshift',shift);
    jcmf=jchmf(par);
    jcmf.psi0=jcmf.MakePsi0*10^-5;
    
end

