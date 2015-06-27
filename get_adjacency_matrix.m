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

