function F_com = centre_of_mass_ln(Z,lattice_dims,n_flux,twist,l,tau)
    %pass the Z = sum(z1,z2,...znp) as Z
    Lx = lattice_dims(1);
    n_particles = size(Z,2);
    q = find_degeneracy(n_particles,n_flux); 
    F_com = thetaln(q*sum(Z,2)/Lx,q*tau,...
        l/q + (n_flux-q)/(2*q)+twist(1)/(2*pi*q),...
        -(n_flux-q)/2-twist(2)/(2*pi));
end