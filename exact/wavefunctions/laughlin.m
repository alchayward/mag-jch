function F = laughlin(Z,lattice_dims,varargin)

if iscell(Z)
    Z = cell2mat(Z).';
end

l  = 0;
if nargin > 2
    l = varargin{1};
end

twist  = [0,0];
if nargin > 3
    twist = varargin{2};
end

q  = 2;
if nargin > 4
    q = varargin{3};
end

Lx = lattice_dims(1);
tau = 1i*lattice_dims(2)/Lx;
n_particles = size(Z,2);
n_flux = q*n_particles;

F = exp(centre_of_mass_ln(Z,lattice_dims,n_flux,twist,l,tau)+...
    landau_level_ln(Z,lattice_dims,n_flux) +...
    laughlin_rel_ln(Z,Lx,tau,q));
end