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

%
% alpha = n_particles/prod(lattice_dims);
% 
% V = perms(1:n_particles);
% laugh = @(xs,l)laughlin(xs,lattice_dims,2,alpha,l,twist,geometry);
% 
% laugh_2 = @(x)laugh(x(1:n_particles/2),l(1)).*...
%                     laugh(x(n_particles/2+1:n_particles),l(2));
% 
% err = 0;
%     function f = pfaff(x)
%         
%         f_all = zeros(1,length(V));
%         for ii = 1:length(V) 
%             [f_all(ii)] = laugh_2(x(V(ii,:)));
%         end
%         f = sum(f_all);
%     end
% 
%     f_out = cellfun(@(x)pfaff(x),X);









