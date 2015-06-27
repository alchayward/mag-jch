function [ psi_list,transition_kappa,max_kappa, diff_list, err_list ] =...
                    psi_vs_kappa( kappa_list,jcmf,jc,varargin )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    opts = struct();
    if nargin > 3
        opts = varargin{1};
    end
    
    if ~isfield(opts,'divirge_kappa')
        opts.divirge_kappa = false;
    end
    divirge_kappa = opts.divirge_kappa;
    
    if ~isfield(opts,'mott_kappa')
        opts.mott_kappa = false;
    end
    mott_kappa = opts.mott_kappa;
    
    psi_list = kappa_list*0;
    diff_list = kappa_list*0;
    err_list = kappa_list*0;
    
    transition_kappa = find_transition_kappa(jc.f_psi,jcmf.amat);
    max_kappa = -jc.onsiteStrength.mu/4;
    for jj = 1:length(kappa_list)
        kappa = kappa_list(jj);
        if kappa > max_kappa && divirge_kappa
            Psi = NaN;
        elseif kappa < transition_kappa && mott_kappa
            Psi = 0;
        else
            [Psi,diff_list(ii),err_list] = jch_lattice_meanfield(...
                kappa,jcmf,jc,opts);
        end
        psi_list(jj) = max(abs(Psi));
        if max(abs(Psi)) > 10^-8
            jcmf.psi0 = Psi; % Use previous value to initate.
        end
    end


