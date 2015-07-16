function F = Pfaffian_L(Z,lattice_dims,varargin)
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

if l == 0
    l1 = 0;
    l2 = 0;
    %twist = twist + [0,pi];
elseif l == 1
    l1 = 1;
    l2 = 1;
    %twist = twist + [pi,pi];
elseif l == 2
    l1 = 1;
    l2 = 0;
    twist = twist + [0,pi];
end

%Lx = lattice_dims(1);
%tau = 1i*lattice_dims(2)/Lx;
n_particles = size(Z,2);
%n_flux = n_particles;
%q=1;

v = perms(1:n_particles);

v1 = v(:,1:n_particles/2);
v2 = v(:,n_particles/2+1:n_particles);
hilbdim = size(Z,1); %it's not really hilb dim, but you know...
F = zeros(hilbdim,1);

n_perms = size(v,1);
for ii = 1:n_perms
    F = F + laughlin(Z(:,v1(ii)),lattice_dims,l1,twist).*...
        laughlin(Z(:,v2(ii)),lattice_dims,l2,twist);
end

end