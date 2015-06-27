lattice_dims = [10;10];
n_points = 100;
butt = zeros(prod(lattice_dims),n_points*2); 
for ii = 0:n_points-1
T = magnetic_lattice_hamiltonian(lattice_dims,.5/n_points*ii);
butt(:,ii+1) = eig(T);
end
butt(:,n_points+1:2*n_points) = fliplr(butt(:,1:n_points));
plot(butt','k.','MarkerSize',5)