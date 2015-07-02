classdef jchmf
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lattice_dims
        n_sites
        amat
        psi0
        alpha
        coords
    end
    
    methods
   
    function obj = jchmf(para)
        obj.lattice_dims = para.lattice_dims;
        obj.alpha = para.alpha(1)/para.alpha(2);
        obj.n_sites=prod(obj.lattice_dims);
        [obj.amat,obj.coords] =...
            magnetic_lattice_hamiltonian(obj.lattice_dims,obj.alpha);
        obj.psi0 = obj.MakePsi0();
    end
        
   
       
   function psi = MakePsi0(obj)
        %[v,~]  = eigs((obj.amat),1,'sr');
        %psi = v(:,1);
        psi = zeros(obj.n_sites,1);
        alpha = obj.alpha;
        for ii =1:obj.n_sites
            X  = obj.coords(:,ii);
            psi(ii) = LatticeSolution(alpha,...
                [1/sqrt(alpha),1/sqrt(alpha)],1i+0.5,X(1)+1i*X(2));
        end
        %psi = rand(length(psi),1);
    end

end
   
    methods(Static)

    function en = energyApprox(psi,a,b,t,amat)
        en = t*abs(psi'*amat*psi)-t^2*a/2*norm(amat*psi)^2+...
            t^4*b/4*norm(amat*psi)^2;
    end

    function psiOut = GSAApprox(psiIn)
        psiOut = -(-1*abs(psiIn).^2+1).*psiIn; 
    end

    end
    
end

