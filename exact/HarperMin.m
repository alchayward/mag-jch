function list = HarperMin( p,q )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


lattice_dims = [p,p];

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

