classdef JCHsystem < Latticesystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

   % energies=struct('kappa',1,'delta',0,'jc',1);
    angle
    atomLevels = 2;
   
        %Operators, like the Jaynes-Cummings interaction, and hopping matrix,
        %are stored as sparse matricies, and only with their upper right
        %triangle elements.

    
        %These lists are Np by Ns matricies, storing the position and number
        %of photonic and atomic occupations at each site.
    end
    
    methods
    
        
        %Constructor
     function obj = JCHsystem(para)   
            obj = obj@Latticesystem(para); 
            if nargin > 0
            end               
            obj.model = 'JCH';
             
     end
          
     
     function p = Initilize(p)
         
         p.hilbDim = p.MakeHilbDim;

         p.siteOccupationList=p.MakeSiteOccupationList;
         
         
         p.sym=p.MakeSymMatrix;
    %     p.angle=p.jchangle;
         p.mat.k=p.MakeKappaMatrix;
         p.mat.d = p.MakeDeltaMatrix;
         p.mat.jc = p.MakeJCMatrix;
         
     end
     
     function h = MakeHamiltonian(obj,varargin)
            nvarargin = length(varargin);
            if nvarargin == 0
                k=obj.hoppingStrength;
                b = obj.onsiteStrength;
            elseif nvarargin == 2
                k = varargin{1};
                b = varargin{2};
            else
                disp(['ERROR - BHystem:MakeHamiltonian:',...
                    'wrong number of arguments']);
            end
            
            h = b(1)*obj.mat.jc + b(2)*obj.mat.d - k*obj.mat.k;
                
            
        end
      
        function h = onsite(obj)
           h = obj.onsiteStrength(1)*obj.mat.jc +... 
               obj.onsiteStrength(3)*obj.mat.d;
        end
        
        function k = hopping(obj)
           k = obj.hoppingStrength(1)*obj.mat.k;
        end

        function d = MakeHilbDim(obj)
            d=0;
            np=obj.nParticles;
            p=partitions(np);
            sites=prod(obj.latticeDim);
            for ii = 1:size(p,1)
                pp=p(ii,:);
                fac=factorial(sum(pp))/prod(factorial(pp));
                s=1;
                for jj=1:np
                   
                    s=s*min([jj+1,obj.atomLevels])^pp(jj);
                end
                d=d+fac*nchoosek(sites,sum(pp))*s;
            end
        end
        
        function list = MakeSiteOccupationList(obj)
         N=obj.nParticles;
         d=obj.nSites;
         dim = obj.hilbDim;     
             %initilize the list to zero.
         global LIST
         LIST = zeros(dim,2*N);
         global II
         II=1;
             
         
             
            function reclist(nn1,dd1,config1)
                if nn1 == 0
                    LIST(II,:) = config1;
                    II = II + 1;
                    II;
                    
                else
                    
                    for jj1 = (dd1):d
                        config1(N-nn1+1) = jj1;
                        reclist2(nn1-1,jj1,config1);
                    end
                    
                    config1(2*N-nn1+1) = 1;
                    
                    
                    for jj1 = dd1:d
                        config1(N-nn1+1) = jj1;
                        reclist2(nn1-1,jj1,config1);
                    end
                   
                end
             end
         

             
        
             function reclist2(nn3,dd3,config3)
                if nn3 == 0
                    LIST(II,:) = config3;
                    II = II + 1;
                    II;
                    
                else
                    
                    config3(N - nn3 + 1) = dd3;
                    reclist2(nn3-1,dd3,config3);
                    
                    for jj3 = dd3+1:d                     
                        config3(N-nn3+1) = jj3;
                        reclist2(nn3-1,jj3,config3)
                    end
                    
                    for jj3 = (dd3+1):d
                        config3(N-nn3+1) = jj3;
                        config3(2*N-nn3+1) = 1;
                        reclist2(nn3-1,jj3,config3)
                    end
                end
             end
             
             
            reclist(N,1,zeros(2*N,1));
         
            list = LIST;
            
         
    
            
         
        end
                     
        function sym = MakeSymMatrix(obj)
  
            
            %for the JCH model, the basis structure is a little more
            %complicated. But we follow a similar proceedure.
            %The basis for distinguishable particles will be 
            %BoxAoxBoxA....Np times. which will be 2^Np times as large as
            %the BH hilbDim.
            
            
            
            
            
            len=factorial(obj.nParticles);
            %posm maps a list of particle positions onto a basis index.
            
            
            %store occuaption as (sss...,aaa...);
            posm=(ones(len,1)*((obj.nParticles-1):-1:0));
            posmb=2*(2*obj.nSites).^posm;
            posma=(2*obj.nSites).^posm;
            
            posm=[posmb';posma']';
            
            
            
            
            list=zeros(obj.hilbDim*len,3);
            iter=1;
            
            %%% The calls to perms is very expensive, so we do it by hand.
            %%% permlist is the permutations to loop over;
            
            permlistb  = perms(1:obj.nParticles);
            permlista  = perms(1+obj.nParticles:2*obj.nParticles);
            permlist = [permlistb,permlista];
            
            
            
            shiftsites = [-ones(1,obj.nParticles),zeros(1,obj.nParticles)];
            
            
            for ii=1:obj.hilbDim
                oclistii = obj.siteOccupationList(ii,:)+shiftsites;
                
                
                v = oclistii(permlist);

              
                
                
              %  li=unique((sum(posm'.*v',1)+1)');
                
                %%%inline unique code
                a = sort(sum(posm'.*v',1)+1);
                d =  logical([ 1 (diff(a)  ~= 0)]);
                li = a(d)';
                
                fac=nnz(li);
                nlist=[ones(1,fac)*ii;li';ones(1,fac)/sqrt(fac)]';
                list(iter:iter+fac-1,:)=nlist;
                iter=iter+fac;
            end
            list=list(1:iter-1,:);
            sym=sparse(list(:,1),list(:,2),list(:,3),...
            obj.hilbDim,(2*obj.nSites)^obj.nParticles);
        end                  
                          
        function d = MakeDeltaMatrix(obj)
               
        a=sparse(kron(eye(obj.nSites),[0,0;0,1]));
        
        d=obj.MakeNParticleMatrix(a);
            
        end 
          
        function jc = MakeJCMatrix(obj)
        a = sparse(kron(eye(obj.nSites),[0,1;1,0]));
        jc=obj.MakeNParticleMatrix(a);
        
        end
                   
        function [tx,ty ] = twist(obj,tangle)
            obj.tangle=tangle;
            [tx,ty] = obj.MakeTwistMatricies;
            
        end

        function li =SingleParticleIndicies(obj,ind)
        
            %for an state space index, return a list of the single particle
            %state indicies.
            
            
            
            occ=obj.siteOccupationList(ind,:);
            pn=obj.nParticles;
            li=zeros(pn,1);
            
            for ii = 1:pn
                li(ii)=1 + (occ(ii)-1)*2 + occ(ii+pn);
            end
        end
            
        function li =SingleParticlePositions(obj,ind)
        
            %for an state space index, return a list of the single particle
            %state indicies.
            
            
            
            occ=obj.siteOccupationList(ind,:);
            pn=obj.nParticles;
            li=zeros(pn,1);
            
            for ii = 1:pn
                pos=obj.sitePositions(occ(ii),:);
                li(ii)=pos(1)+1i*pos(2);
                
            end
        end
          
        function [mx,my,mtx,mty] = MakeAdjacencyMatricies(obj)
        
            b=[1,0;0,0];
             switch obj.topology
                case 'torus'
                    [ mx,my ] = obj.MakeHardAdjacencyMatrices;
                    [ mtx,mty ] = obj.MakeWrapAdjacencyMatrices;
                    mx=kron(mx,b);
                    my=kron(my,b);
                    mtx=kron(mtx,b);
                    mty=kron(mty,b);
                case 'hard'
                    [ mx,my ] = obj.MakeHardAdjacencyMatrices;
                    mtx = [];
                    mty = [];
                    my=kron(my,b);
                    mx=kron(mx,b);
             end
                    
        end
        
        function  f = posfac(obj,ii)
                  %for each particle, if in boson/atomic state, factor is
                  %given by cos/sin(th(alpha,beta))
                 
                  np=obj.nParticles;
                  fac=obj.angle;
                  p=obj.siteOccupationList(ii,(np+1):(2*np));
                  f=1;
                  
                  for ii = p 
                    f=f*fac(ii+1);
                  end
                  
                  
        end   
        
        function fac = jchangle(obj,beta,delta,kappa)
       
                  
                  if obj.alpha ==1 || obj.alpha == 0
                      k = -4;
                  else
                    n=round(abs(obj.alpha)*prod(obj.latticeDim));
                    e=HarperMin(prod(obj.latticeDim));
                    k=e(n,2);
                  end
                  ssite=[0,beta;...
                      beta,delta-...
                      kappa*k];
                  
                  [v,~]=eig(ssite); 
                  fac=v(:,1);
        end
        
        function pro = PNProjector(obj,pn)
           p=obj.nParticles;
           pro=zeros(1,obj.hilbDim);
           if pn ==1
            for ii = 1:obj.hilbDim
               if length(unique(nonzeros(obj.siteOccupationList(ii,1:p))))==p
                   pro(ii)=1;
               end
            end
             pro=diag(sparse(pro));
           elseif pn==2
               pro=sparse(1:obj.hilbDim,1:obj.hilbDim,1)-obj.PNProjector(1);
           elseif pn==3
           for ii = 1:obj.hilbDim
               if length(unique(nonzeros(obj.siteOccupationList(ii,1:p))))==1
                   pro(ii)=1;
                   
               end
           end
           pro=diag(sparse(pro));
           end
           
           
        end
        
        function e = E0est(obj)
        
         [p,q] = rat(obj.alpha);
         hop=obj.hoppingstrength*HarperMin(q,p);
         e=obj.nParticles*jchhe(1,1,obj.onsitestrength(2),hop,-1);
        end
    
        function factors = atomic_state_factor(obj,dims,varargin)
            beta = obj.onsiteStrength(1);
            delta = obj.onsiteStrength(2);
            kappa = obj.onsiteStrength(1);
            if nargin > 2
                params = varargin{1};
                if isfield(params,'beta')
                    beta = params.beta;
                end
                if isfield(params,'delta')
                    delta = params.delta;
                end
                if isfield(params,'kappa')
                    kappa = params.kappa;
                end
            end
            np =obj.nParticles;
            ang = obj.jchangle(beta,delta,kappa);
            %ss is the total number of atomic excitations for each dim.
            ss =sum(obj.siteOccupationList(dims,np+1:2*np),2);
            %a factor to multiply each state by. 
            factors = (ang(2).^ss).*(ang(1).^(np-ss));
        end
        
        function psiL = LaughlinState(obj,varargin)
            
            params = struct();
            twist = obj.tangle;
            if nargin > 1
                params = varargin{1};
                if isfield(params,'twist')
                    twist = params.twist;
                end
            end
            facs = obj.atomic_state_factor(1:obj.hilbDim,params);
            np = obj.nParticles;
            Z = reshape(...
                obj.sitePositions(obj.siteOccupationList(:,1:np)',1)+...
                obj.sitePositions(obj.siteOccupationList(:,1:np)',2)*1i,...
                np,obj.hilbDim).';
            
            lattice_dims = obj.latticeDim;
            
            psiL = zeros(obj.hilbDim,2);
            wf = Wavefunctions('laughlin');
            psiL(:,1) = facs.*wf(Z,lattice_dims,1,twist);
            psiL(:,2) = facs.*wf(Z,lattice_dims,2,twist);

            h = HubbardLibrary();
            psiL = h.GramSchmit(psiL);
            
        end
    end %methods
    
end %classdef

