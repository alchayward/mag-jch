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