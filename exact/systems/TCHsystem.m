classdef TCHsystem < Latticesystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

   % energies=struct('kappa',1,'delta',0,'jc',1);
    angle
    atomLevels = [2 1 1];
    nAtoms = 2;
   
        %Operators, like the Jaynes-Cummings interaction, and hopping matrix,
        %are stored as sparse matricies, and only with their upper right
        %triangle elements.

    
        %These lists are Np by Ns matricies, storing the position and number
        %of photonic and atomic occupations at each site.

    
    end
    
    methods
    
        
        %Constructor
     function obj = TCHsystem(para)
            
            obj = obj@Latticesystem(para); 
            if nargin > 0
                   
            end               
            obj.model = 'TCH';
            obj.nAtoms = para.nAtoms;

	    %atomLevels(1) = obj.nParticles;
	    %obj.atomLevels = para.atomLevels;
        obj.atomLevels(1) = para.maxParticlesPerSite;
     end
          
     
     function p = Initilize(p)
         
         p.hilbDim = 1;

         p.siteOccupationList=p.MakeSiteOccupationList;
         p.hilbDim = length(p.siteOccupationList);
         
         p.sym=p.MakeSymMatrix;
         p.mat.k=p.MakeKappaMatrix;
         [p.mat.d1 p.mat.d2] = p.MakeDeltaMatrix;
         [p.mat.jc1 p.mat.jc2 p.mat.jc3] = p.MakeJCMatrix;
         
     end
     
     
     
     
     
     
      %System Things
      
    
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
            
            h = b(1)*obj.mat.jc1 + b(2)*obj.mat.jc2 + b(3)*obj.mat.jc3 ...
                + b(4)*obj.mat.d1 + b(5)*obj.mat.d2...
                - k*obj.mat.k;
                
            
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
         
            function reclist(nn,dd,config,atLimit)
                if nn == 0
                    LIST(II,:) = config;
                    II = II + 1;
                    II;
                    
                else
                    for ii = 1:(obj.nAtoms + 1)
                        
                        config(2*N-nn+1) = ii-1;
                        atLimitTemp = atLimit;
                        atLimitTemp((ii+1):(obj.nAtoms+1))=0;
                        atLimitTemp(ii)=atLimitTemp(ii)-1;
                      
                        if atLimit(ii) > 0
                            config(N-nn+1) = dd;
                            %atLimitTemp = atLimit;
                           % atLimitTemp((ii+1):(obj.nAtoms+1))=0;
                          %  atLimitTemp(ii)=atLimitTemp(ii)-1;
                            reclist(nn-1,dd,config,atLimitTemp);
                        end
                        
                       % atLimitTemp = atLimit;
                       % atLimitTemp(ii)=atLimitTemp(ii)-1;
                        
                        for jj = (dd + 1):d
                            config(N-nn+1) = jj;
                            reclist(nn-1,jj,config,atLimitTemp);
                        end
                    end
                        
                end
            end
         
            
              
             
            reclist(N,1,zeros(2*N,1),obj.atomLevels);
         
            list = LIST;
         
        end
                     
        function sym = MakeSymMatrix(obj)
  
            
            %for the JCH model, the basis structure is a little more
            %complicated. But we follow a similar proceedure.
            %The basis for distinguishable particles will be 
            %BoxAoxBoxA....Np times. which will be Na^Np times as large as
            %the BH hilbDim.
            
            
            np = obj.nParticles;
            ns = prod(obj.latticeDim);
            hd = obj.hilbDim;
            na = obj.nAtoms+1;
            symDim = ((na)*ns)^np;
            
            len=factorial(np);
            %posm maps a list of particle positions onto a basis index.
            
            
            %store occuaption as (sss...,aaa...);
            posm=(ones(len,1)*((np-1):-1:0));
            posmb=na*(na*ns).^posm;
            posma=(na*ns).^posm;
            
            posm=[posmb';posma']';

            list=zeros(hd*len,3);
            iter=1;
            
            %%% The calls to perms is very expensive, so we do it by hand.
            %%% permlist is the permutations to loop over;
            
            permlistb  = perms(1:np);
            permlista  = perms(1+np:2*np);
            permlist = [permlistb,permlista];
            %nperms=size(permlist,1);
            
            %v = permlist;
            shiftsites = [-ones(1,np),zeros(1,np)];
            
            
            for ii=1:hd
                oclistii = obj.siteOccupationList(ii,:)+shiftsites;
                
                
                v = oclistii(permlist);

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
            hd,symDim);
        end                  
                          
        function [d1 d2] = MakeDeltaMatrix(obj)
               
        d1o = [0,0,0;0,1,0;0,0,0];
        d2o = [0,0,0;0,0,0;0,0,1];
        a=sparse(kron(eye(obj.nSites),d1o));
        
        d1=obj.MakeNParticleMatrix(a);
        
        a=sparse(kron(eye(obj.nSites),d2o));
        
        d2=obj.MakeNParticleMatrix(a);
            
        end 
          
        function [jc1 jc2 jc3] = MakeJCMatrix(obj)
        
        jc1o = [0,1,0;1,0,0;0,0,0];
        jc2o = [0,0,1;0,0,0;1,0,0];
        jc3o = [0,0,0;0,0,1;0,1,0];
        
        a=sparse(kron(eye(obj.nSites),jc1o));
        jc1=obj.MakeNParticleMatrix(a);
         
        a=sparse(kron(eye(obj.nSites),jc2o));
        jc2=obj.MakeNParticleMatrix(a);
        
        a=sparse(kron(eye(obj.nSites),jc3o));
        jc3=obj.MakeNParticleMatrix(a);
        
        end
                   
        function [tx ty ] = twist(obj,tangle)
            obj.tangle=tangle;
            [tx ty] = obj.MakeTwistMatricies;
            
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
          
        function [mx my mtx mty] = MakeAdjacencyMatricies(obj)
        
            b=[1,0,0;0,0,0;0,0,0];
             switch obj.topology
                case 'torus'
                    [ mx my ] = obj.MakeHardAdjacencyMatrices;
                    [ mtx mty ] = obj.MakeWrapAdjacencyMatrices;
                    mx=kron(mx,b);
                    my=kron(my,b);
                    mtx=kron(mtx,b);
                    mty=kron(mty,b);
                case 'hard'
                    [ mx my ] = obj.MakeHardAdjacencyMatrices;
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
        
        function fac = jchangle(obj)
       
                  n=round(abs(obj.alpha)*prod(obj.latticeDim));
                  e=HarperMin(prod(obj.latticeDim));
                  if obj.alpha ==1 || obj.alpha == 0
                      k = -4;
                  else
                    k=e(n,2);
                  end
                  ssite=[0,obj.onsitestrength(1);...
                      obj.onsitestrength(1),obj.onsitestrength(2)-...
                      obj.hoppingstrength*k];
                  
                  [v ~]=eig(ssite); 
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
        
         [p q] = rat(obj.alpha);
         hop=obj.hoppingstrength*HarperMin(q,p);
         
  
         e=obj.nParticles*jchhe(1,1,obj.onsitestrength(2),hop,-1);
     end
    
        
        
    end %methods
    
   
    
    
end %classdef

