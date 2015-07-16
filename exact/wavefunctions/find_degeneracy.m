function q = find_degeneracy(n_particles,n_flux)

g = gcd(n_particles,n_flux);
q = n_flux/g;

end