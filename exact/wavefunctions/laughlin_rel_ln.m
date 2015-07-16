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