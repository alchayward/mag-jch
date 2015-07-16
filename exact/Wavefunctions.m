function h = Wavefunctions()
%HubbardLibrary: Library of functions for Hubbard Model Calculation
%   Functions: Chern(systemParameters, chernParameters)
%              
h.Pfaffian = @Pfaffian;
h.Laughlin = @Laughlin;
h.landau_ln = @landau_ln;
h.pfaffian_rel_ln = @pfaffian_rel_ln;
h.pfaffian_rel_ln = @pfaffian_rel_ln;
h.pfaffian_matricies = @pfaffian_matricies;
h.perm_sign = @perm_sign;
h.recursive_pairings = @recursive_pairings;
h.pair_ind = @pair_ind;
h.laughlin_rel_ln = @laughlin_rel_ln;
h.centre_of_mass_ln = @centre_of_mass_ln;
h.find_degeneracy = @find_degeneracy;
h.find_magnetic_length = @find_magnetic_length;
h.landau_level_ln = @landau_level_ln;
h.Laughlin = @Laughlin;
end

function F = pfaffian_rel_ln(Z,Lx,tau,l,q)
n_particles = size(Z,2);

th_1_func = @(z)thetaln(z/Lx,tau,1);
if l == 0
    %theta_2
    th_a_func = @(z)thetaln(z/Lx,tau,2);
elseif l == 1
    %theta_3
    th_a_func = @(z)thetaln(z/Lx,tau,3);
elseif l == 2
    %theta_4
    th_a_func = @(z)thetaln(z/Lx,tau,4);
end


pairings = laughlin_pairs(n_particles);
n_pairings = size(pairings,1);
[pfaffian_pairs, pfaffian_signs] = pfaffian_matricies(n_particles);
n_pfaff_pairings = size(pfaffian_pairs,1);

Z_diffs = Z(:,pairings(:,1))-Z(:,pairings(:,2));

hilb_dim = size(Z,1);
F_rel = zeros(hilb_dim,1);

for ii = 1:n_pfaff_pairings %sum over each pfaffian pairing    
    pairs = reshape(pfaffian_pairs(ii,:).',2,n_particles/2).';
    pfaff_inds = zeros(n_particles/2,1);
    for p=1:n_particles/2
        pfaff_inds(p) = pair_ind(pairs(p,:),pairings);
    end
    
    t1_inds = 1:n_pairings;
    for jj = pfaff_inds.'
        t1_inds = t1_inds(t1_inds~=jj);
    end
    F_rel_ln = zeros(hilb_dim,1);
    
    for jj =1:n_particles/2
        F_rel_ln = F_rel_ln + th_a_func(Z_diffs(:,pfaff_inds(jj)));
    end
       
    for jj =1:(n_pairings-n_particles/2)
        F_rel_ln = F_rel_ln + q*th_1_func(Z_diffs(:,t1_inds(jj)));
    end
    F_rel = F_rel + pfaffian_signs(ii)*exp(F_rel_ln);
end


F = log(F_rel);

end

function [pairs, signs] = pfaffian_matricies(n_particles)

%recursivly generate all pairings
        pairs = recursive_pairings(1:n_particles);
        signs = zeros(size(pairs,1),1);
        for ii = 1:size(pairs,1)
            signs(ii) = perm_sign(pairs(ii,:));
        end
end
    
function sig = perm_sign(p)
I = speye(length(p));
sig = det(I(:,p));           
end

function pairs_out = recursive_pairings(particles)
np = length(particles);
if np == 2
        pairs_out = particles;
else
    p1 = particles(1);
    particles = particles(particles~=p1);
    pairs_out = zeros(factd(np-1),np);
    ind = 1;
    for p2 = particles

        d_particles = particles(particles~=p2);
        pair = [p1,p2];
        x = recursive_pairings(d_particles);
        
        for jj = 1:size(x,1);
            pairs_out(ind,:) = [pair,x(jj,:)];
            ind = ind+1;
        end
    end
end

end

function ind = pair_ind(pair,pair_list)
   ind = 0;
    for pl = 1:size(pair_list,1)
       if isequal(pair,pair_list(pl,:))
           ind = pl;
       end
    end
end

function F_rel_ln = laughlin_rel_ln(Z,Lx,tau,q)
n_particles = size(Z,2);
pairings = laughlin_pairs(n_particles);
n_pairings = size(pairings,1);
Z_diffs = Z(:,pairings(:,1))-Z(:,pairings(:,2));
hilb_dim = size(Z,1);

F_rel_ln = zeros(hilb_dim,1);
for ii = 1:n_pairings
    F_rel_ln = F_rel_ln + q*thetaln(Z_diffs(:,ii)/Lx,tau,1);
end

function pairings = laughlin_pairs(n_particles)
n_pairings = n_particles*(n_particles-1)/2;
pairings = zeros(n_pairings,2);
ind = 1;
for ii = 1:n_particles-1
    for jj = ii+1:n_particles
        pairings(ind,:) = [ii,jj];
        ind = ind + 1;
    end
end

function F = Laughlin(Z,lattice_dims,varargin)

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

function F_landau_ln = landau_level_ln(Z,lattice_dims,n_flux)
    Lb2 = find_magnetic_length(lattice_dims,n_flux)^2;
    F_landau_ln = -sum(imag(Z).^2,2)/(2*Lb2);
end

function Lb = find_magnetic_length(lattice_dims,n_flux)
    Lb = sqrt(prod(lattice_dims)/(2*pi*n_flux));
end

function q = find_degeneracy(n_particles,n_flux)
g = gcd(n_particles,n_flux);
q = n_flux/g;
end

function F_com = centre_of_mass_ln(Z,lattice_dims,n_flux,twist,l,tau)
    %pass the Z = sum(z1,z2,...znp) as Z
    Lx = lattice_dims(1);
    n_particles = size(Z,2);
    q = find_degeneracy(n_particles,n_flux); 
    F_com = thetaln(q*sum(Z,2)/Lx,q*tau,...
        l/q + (n_flux-q)/(2*q)+twist(1)/(2*pi*q),...
        -(n_flux-q)/2-twist(2)/(2*pi));
end

function F = Pfaffian(Z,lattice_dims,varargin)
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
    %theta_2
    twist = twist + [0,pi];
elseif l == 1
    %theta_3
    twist = twist + [pi,pi];
elseif l == 2
    %theta_4
    twist = twist + [pi,0];
end

Lx = lattice_dims(1);
tau = 1i*lattice_dims(2)/Lx;
n_particles = size(Z,2);
n_flux = n_particles;
q=1;

F = exp(centre_of_mass_ln(Z,lattice_dims,n_flux,twist,l,tau)+...
    landau_level_ln(Z,lattice_dims,n_flux) +...
    pfaffian_rel_ln(Z,Lx,tau,l,q));
end
