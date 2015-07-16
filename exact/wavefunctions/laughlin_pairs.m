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