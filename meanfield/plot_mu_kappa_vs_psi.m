function [ h ] = plot_mu_kappa_vs_psi(M,mu_list,kappa_list,...
    varargin)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

plot_mott = false;
if nargin > 3
    trans_list = varargin{1};
    plot_mott = true;
end

plot_divirge = false;
if nargin > 4
    divirge_list = varargin{2};
    plot_divirge = true;
end



h = imagesc(kappa_list,mu_list,flipud(M));

hold on
colorbar
%%%Mott region
if plot_mott
    mott_M = (1./trans_list)'*kappa_list<1;
    green = cat(3, zeros(size(mott_M)),ones(size(mott_M)), zeros(size(mott_M)));
    hg = imagesc(kappa_list,mu_list,green); 
    set(hg, 'AlphaData', flipud(mott_M));
    plot(trans_list,fliplr(mu_list),'k');
end
%%%Divirgent region
if plot_divirge
    divirge_M = (1./divirge_list)'*kappa_list>1;
    white = cat(3, ones(size(mott_M)),ones(size(mott_M)), ones(size(mott_M)));
    hg = imagesc(kappa_list,mu_list,white); 
    set(hg, 'AlphaData', flipud(divirge_M));
    plot(divirge_list,fliplr(mu_list),'k');
end
hold off
end

