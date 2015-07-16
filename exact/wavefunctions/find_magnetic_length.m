function Lb = find_magnetic_length(lattice_dims,n_flux)
    Lb = sqrt(prod(lattice_dims)/(2*pi*n_flux));
end