function [ T ] = magnetic_lattice_hamiltonian(...
        lattice_dims, alpha, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%This only works when A_func is linear!!!!

lattice_dims = reshape(lattice_dims,2,1);
shift_vec = [0;0];
X_wrap = @(X)[mod(X(1),lattice_dims(1));mod(X(2),lattice_dims(2))] + ...
                shift_vec;
dim = prod(lattice_dims);

default_vector_potential = 'landau';
A_func = get_vector_potential_func(default_vector_potential,alpha);
if nargin > 3
    if isa(varargin{2},'char')
        A_func = get_vector_potential_func(varargin{2},varargin{3});
    elseif isa(varargin{2},'function_handle')
        A_func = varargin{2};
    else
        disp('invalid vector potential, defaulting');
    end
end

PBC = [0;0];
boundary_conditions = 'torus';
if nargin > 3
    %check that PCB is a length 2 vector
    if isequal(size(varargin{1}),[2,1]);
        PBC = varargin{1};
    else
        disp('invalid boundary conditions');
    end
end

is_torus = false;
if strcmp(boundary_conditions,'torus')
    is_torus = true;
end

Mag_BC = @(X)(2*pi*alpha*[X(2);X(1)]-A_func(X)).*lattice_dims;


%lattice positions
n_sites = prod(lattice_dims);

[X_mat,Y_mat] = meshgrid(0:lattice_dims(1)-1,0:lattice_dims(2)-1);
coords = [reshape(X_mat,n_sites,1),reshape(Y_mat,n_sites,1)]';

%invert the coords (ie, find the index corresponding to a position)
%These 1s due to the array indexing are pretty annoying.
coords_ind = zeros(lattice_dims(1),lattice_dims(2));
for ii = 1:dim
    coords_ind(coords(1,ii)+1,coords(2,ii)+1) = ii;
end

udlr = [0,1;0,-1;-1,0;1,0]';
T = zeros(dim);
for ii = 1:dim
    X = coords(:,ii);
    for dx = udlr
        Y = X+dx;
        wrap = any(floor(Y./lattice_dims));
        Y = X_wrap(Y);
        jj = coords_ind(Y(1)+1,Y(2)+1);
        phase = dx.'*(A_func(X));
        if wrap && is_torus
            bc = dx'*(PBC-Mag_BC(X));
            phase = phase + bc;
        end
        T(ii,jj) = exp(-1i*phase);
    end
end

end

