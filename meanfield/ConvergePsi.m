function [ psi0 max_diff, err] = ConvergePsi(psi0, f_Psi, varargin)
    %ConvergePsi: take in a vector of order parameters, and a function
    %f_Psi(psi), which will return a set of order parameters.
    % opts is a struct with possible fields(default):
        % max_iter(1000) : maximum interations
        % min_error(10^-10) : convergence criteria
        % iter_plot (False) : plot something each iteration
        % 

    %Set some default paramters for the funcion if not provided
    opts = struct();
    if nargin > 2
        opts = varargin{1};
    end

    if ~isfield(opts,'max_iter')
        opts.max_iter = 1000;
    end
    max_iter = opts.max_iter;
    
    if ~isfield(opts,'min_diff')
        opts.min_diff = 10^-10;
    end
    min_diff = opts.min_diff;
    
    if ~isfield(opts,'mod_psi')
        opts.mod_psi = 1;
    end
    mod = opts.mod_psi;
    
    graph_it = false;
    lattice_dims = [];
    if isfield(opts,'iter_plot')
        if isfield(opts,'lattice_dims')
           if length(opts.lattice_dims) == 2
               graph_it = opts.iter_plot;
               lattice_dims = opts.lattice_dims;
           end
        end
    end

   if graph_it
     %  hline=plot(abs(psi0));
      hline=imagesc(abs(reshape(psi0, lattice_dims(1),...
          lattice_dims(2))));
      shading flat;
      colorbar
   end
   %title(num2str(obj.onsiteStrength.k))

   ii=0;
   err = 0;
   max_diff = min_diff+1;
   while ((ii < max_iter) && (abs(max_diff) > min_diff))
        a = f_Psi(psi0);
        psi = psi0*(1-mod) + mod*a; 
        psiDiff = psi-psi0; 
        max_diff = max(abs(psiDiff));%/max(abs(psi));

        if graph_it
            set(hline,'CData',abs(reshape(psi,lattice_dims(1),...
                lattice_dims(2))));  %# Update the y data of the line
            xlabel(['diff = ' , num2str(max_diff)])
            ylabel(num2str(ii))
            drawnow       %# Force the graphics to update immediately
        end
        psi0=psi;
        ii=ii+1;
        if max_diff == inf
            err = 2;
            fprintf('iteration divirged: iter_num = %d\n',ii);
            break
        end
   end
    
    if ii >= max_iter
         fprintf('failed to converge!: diff = %s\n',max_diff);
         err = 1;    
    end
    
end

