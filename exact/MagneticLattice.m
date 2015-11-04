function [ ml ] = MagneticLattice()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ml = struct();
ml.magnetic_lattice_hamiltonian_XY = @magnetic_lattice_hamiltonian_XY;
ml.magnetic_lattice_hamiltonian = @magnetic_lattice_hamiltonian;
ml.get_vector_potential = @get_vector_potential;
ml.get_adjacency_matrix = @get_adjacency_matrix;
ml.HarperMin = @HarperMin;
ml.thetaln = @thetaln;
ml.theta3ln = @theta3ln;
ml.theta = @theta;

end

function [ TX, TY ] = magnetic_lattice_hamiltonian_XY(...
        lattice_dims, alpha, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%This only works when A_func is linear!!!!

%%%This should be cominded with magnetic_lattice_hamiltonian and 
%%%Diffirenitated through inouts and outputs.

lattice_dims = reshape(lattice_dims,2,1);
shift_vec = [0;0];
X_wrap = @(X)[mod(X(1),lattice_dims(1));mod(X(2),lattice_dims(2))] + ...
                shift_vec;
dim = prod(lattice_dims);
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

default_vector_potential = 'landau';
A_func = get_vector_potential(default_vector_potential,alpha);
if nargin > 3
    if isa(varargin{2},'char')
        A_func = get_vector_potential(varargin{2},varargin{3});
    elseif isa(varargin{2},'function_handle')
        A_func = varargin{2};
    else
        disp('invalid vector potential, defaulting');
    end
end

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


%adjaceny matrix (1s and 0s)
%%% This is probably inefficient, since we're iterating over all the 
%%% edges twice.
[a_mat,b_mat] = get_adjacency_matrix(coords,boundary_conditions);
a_mat = logical(a_mat);
b_mat = logical(b_mat);

udlr = [0,1;0,-1;-1,0;1,0]';
TX = zeros(dim);
TY = TX;

%Mag_BC = @(X)(2*pi*alpha*[X(2);X(1)]-A_func(X)).*lattice_dims;
Mag_BC = @(X)(2*pi*alpha*[X(2);X(1)]-A_func(X)).*lattice_dims;


for ii = 1:dim
    X = coords(:,ii);
    dx = udlr(:,1);
    Y = X_wrap(X+dx);
    jj = coords_ind(Y(1)+1,Y(2)+1);
    phase = dx.'*(A_func(X));
    if b_mat(ii,jj)
        bc = dx'*(PBC-Mag_BC(X));
        phase = phase + bc;
    end
    TY(ii,jj) = exp(-1i*phase);
    
    dx = udlr(:,4);
    Y = X_wrap(X+dx);
    jj = coords_ind(Y(1)+1,Y(2)+1);
    phase = dx.'*(A_func(X));
    if b_mat(ii,jj)
        bc = dx'*(PBC-Mag_BC(X));
        phase = phase + bc;
    end
    TX(ii,jj) = exp(-1i*phase);
    
end
TX = TX.*(a_mat+b_mat);
TY = TY.*(a_mat+b_mat);
end

function [ T, coords ] = magnetic_lattice_hamiltonian(...
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
A_func = get_vector_potential(default_vector_potential,alpha);
if nargin > 3
    if isa(varargin{2},'char')
        A_func = get_vector_potential(varargin{2},varargin{3});
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
            boundary_conditions = dx'*(PBC-Mag_BC(X));
            phase = phase + boundary_conditions;
        end
        T(ii,jj) = exp(-1i*phase);
    end
end

end

function [ A_func ] = get_vector_potential( gauge, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if strcmp(gauge,'landau')
    A_func = @(X)2*pi*alpha*[0;X(1)];
elseif strcmp(gauge,'symmetric')
    A_func = @(X)2*pi*alpha*[-X(2);X(1)]/2;
end
end

function [ a_mat,b_mat ] = get_adjacency_matrix(coords,...
                                        boundary_conditions )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

X = coords(1,:);
Y = coords(2,:);
dim = length(X);
a_mat = zeros(dim);
sep_bound = false;
if nargout == 2
    sep_bound = true;
    b_mat = a_mat;
end

isboundary = @(M,ii,jj)M(ii) == min(M) && M(jj) == max(M); 

    for ii = 1:dim
        for jj = 1:dim
            diff_x = abs(X(ii)-X(jj));
            diff_y = abs(Y(ii)-Y(jj));
            if (diff_x == 1 && diff_y == 0) || (diff_x == 0 && diff_y == 1)
                a_mat(ii,jj) = 1;
            end
        
            if strcmp(boundary_conditions,'torus')
                if (isboundary(X,ii,jj) || isboundary(X,jj,ii)) && Y(ii)==Y(jj)
                    if sep_bound
                        b_mat(ii,jj) = 1;
                    else
                        a_mat(ii,jj) = 1;
                    end
                elseif (isboundary(Y,ii,jj) || isboundary(Y,jj,ii)) && X(ii)==X(jj)
                    if sep_bound
                        b_mat(ii,jj) = 1;
                    else
                        a_mat(ii,jj) = 1;
                    end
                end
            end
        end
    end
end

function list = HarperMin( p,q )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
lattice_dims = [p,3];

if nargin == 2
    alpha=(q)/p;
    A=magnetic_lattice_hamiltonian(lattice_dims, alpha);
    list=min(eig(full(A)));
else
    list=zeros(p-1,2);
    for ii=1:p-1
        alpha=(ii)/p;
        A=magnetic_lattice_hamiltonian(lattice_dims, alpha);
        list(ii,:)=[ii/p,min(eig(full(A)))];
    end
end
end

function f = thetaln(z,tau,varargin)
% Maps the genrealized Jacobi elliptical theta function to theta3. Returns
% the log of theta[a,b](z|tau)
iter_lim = 7;

if nargin == 3 %take in a standard elliptic theta (1,2,3,4)
    switch varargin{1}
        case 1
            %theta_1
            a = 1/2;
            b = 1/2;
        case 2
            %theta_2
            a = 1/2;
            b = 0;
        case 3
            %theta_3
            a = 0;
            b = 0;
        case 4
            %theta_4
            a = 0;
            b = 1/2;
    end
elseif nargin == 4
    a = varargin{1};
    b = varargin{2};
end
    
a0 = mod(a,1);
z1=z+b+a0*tau;
f_mod = 1i*pi*(tau*a0^2+2*a0*(z+b));

%remove periodic elements so that z lies in the (0,0) to (1,tau) reigion.
%This makes gaurentees about convergence better. Assuming tau = 1 means
%iter_lim = 7 is very well converged.
n = floor(imag(z1)./imag(tau));
z2 = z1 - n*tau;
m = floor(real(z2));
z3 = z2 - m;

iter = -iter_lim:iter_lim;
f_sum = sum(exp(1i*pi*tau*ones(length(z),1)*(iter.^2)+2i*pi*z3*iter),2);

f = -1i*pi*(tau*(n.^2)+2*n.*z3)+log(f_sum)+f_mod;
end

function f = theta3ln(z,tau)
iter_lim = 7;

%remove periodic elements so that z lies in the (0,0) to (1,tau) reigion.
%This makes gaurentees about convergence better. Assuming tau = 1 means
%iter_lim = 7 is very well converged.
n = floor(imag(z)./imag(tau));
z1 = z - n*tau;
m = floor(real(z1));
z2 = z1 - m;

iter = -iter_lim:iter_lim;
f_sum = sum(exp(1i*pi*tau*ones(length(z),1)*(iter.^2)+2i*pi*z2*iter),2);

f = -1i*pi*(tau*(n.^2)+2*n.*z2)+log(f_sum);
end

function f = theta(a,b,z,tau)

iterMax = 200;
inds = -iterMax:iterMax;
theta_fun = @(n)exp(1i*pi*tau*(a+n)^2+2i*pi*(n+a)*(z+b));
f = sum(arrayfun(theta_fun,inds));
end