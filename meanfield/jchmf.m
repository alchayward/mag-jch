classdef jchmf < Latticesystem
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        dim
        amat
        psi0
    end
    
    methods
   
    function obj = jchmf(para)
        para.nParticles=1;
        para.hoppingStrength=1;
        para.A = [-1,0];
        para.maxParticlesPerSite=1;
        para.sysPath='./vorticies/';
        para.twist = 0*[pi/2,pi/6];

        obj = obj@Latticesystem(para); 

        obj.nSites=obj.NumberOfSites;
        obj.sitePositions=obj.MakeSitePositions;
        obj.amat = obj.MakeAdjacencyMatrix;
    end
        
   
       
   function psi = MakePsi0(obj)
        %[v,~]  = eigs((obj.amat),1,'sr');
        %psi = v(:,1);
        psi = zeros(1,obj.nSites);
        alpha = obj.alpha;
        for ii =1:obj.nSites
            X  = obj.sitePositions(ii,:);
            psi(ii) = LatticeSolution(alpha,...
                [1/sqrt(alpha),1/sqrt(alpha)],1i+0.5,X(1)+1i*X(2));
        end
        psi = psi';
        %psi = rand(length(psi),1);
    end

   function tm = MakeTrapMat(obj)
        centre = (obj.latticeDim + 1)/2;
        tm = diag(sum((obj.sitePositions-ones(obj.nSites,1)*centre).^2,2));
   end

   function a = MakeAdjacencyMatrix(obj)
       %[mx, my] = MakeHardAdjacencyMatrices(obj);
       %[mtx, mty] = MakeWrapAdjacencyMatrices(obj);
       %a = -(mx + my + mtx + mty);
       %a=a+a';
       %alpha_div = alpha(1)/alpha(2);
       a = magnetic_lattice_hamiltonian(obj.latticeDim,obj.alpha);
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

