function F_landau_ln = landau_level_ln(Z,lattice_dims,n_flux)
    Lb2 = find_magnetic_length(lattice_dims,n_flux)^2;
    F_landau_ln = -sum(imag(Z).^2,2)/(2*Lb2);
end