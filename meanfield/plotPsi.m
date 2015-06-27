function hline = plotPsi(x, dim)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    hline=imagesc(abs(reshape(x, dim(1),...
        dim(2))));
    shading flat;
    colorbar
end

