classdef GJCHsystem < Latticesystem
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
     function obj = GJCHsystem(para)
            
            obj = obj@Latticesystem(para); 
            if nargin > 0
            obj.atomLevels = para.atomLevels;       
            end               
            obj.model = 'GJCH';
            
             
     end
          
     
     function p = Initilize(p)
         
         p.hilbDim = p.MakeHilbDim;

         p.siteOccupationList=p.MakeSiteOccupationList;
         
         
         p.sym=p.MakeSymMatrix;
    %     p.angle=p.jchangle;
         p.mat.k=p.MakeKappaMatrix;
         
          
         p.mat.d = {};
         for ii = 1:(p.atomLevels-1)
            p.mat.d{ii} = p.MakeDeltaMatrix(ii);
         end
         
         
         
         p.mat.jc = {};
         for ii = 1:(p.atomLevels-1)
            p.mat.jc{ii} = p.MakeJCMatrix(ii);
         end
         
         
        % [p.mat.d1 p.mat.d2] = p.MakeDeltaMatrix;
        % [p.mat.jc1 p.mat.jc2] = p.MakeJCMatrix;
         
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
            
            h= - k*obj.mat.k;
            for ii = 1:(obj.atomLevels-1)
                h= h + b(ii)*obj.mat.jc{ii} +...
                    b(ii+obj.atomLevels-1)*obj.mat.d{ii};
            end
            
            
        end
      
      
      
        function h = onsite(obj)
           h = obj.onsitestrength(1)*obj.mat.jc1 +... 
               obj.onsitestrength(2)*obj.mat.jc2 +...
               obj.onsitestrength(3)*obj.mat.d1 + ...
               obj.onsitestrength(4)*obj.mat.d2;
        end
        
        function k = hopping(obj)
           k = obj.hoppingstrength(1)*obj.mat.k;
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
                    if pp(jj) < obj.atomLevels
                        
                    end
                    s=s*min([jj+1,obj.atomLevels])^pp(jj);
                end


                d=d+fac*nchoosek(sites,sum(pp))*s;

            end

        end
        
        function list = MakeSiteOccupationList(obj)

        
         N=obj.nParticles;
         d=obj.nSites;
         dim = obj.hilbDim;
         at0=obj.atomLevels-1;
         
        
     
     
           
         
             
             %initilize the list to zero.
         global LIST
         LIST = zeros(dim,2*N);
         global II
         II=1;
         
         
           
         
             
            function reclist(nn1,dd1,config1,at)
                if nn1 == 0
                    LIST(II,:) = config1;
                    II = II + 1;
                    II;
                    
                else
                   

                    for jj1 = (dd1):d
                        config1(N-nn1+1) = jj1;
                        reclist(nn1-1,jj1,config1,0);
                    end

                    config1(2*N-nn1+1) = 1;

                     for jj1 = (dd1+1):d
                        config1(N-nn1+1) = jj1;
                        reclist(nn1-1,jj1,config1,at0-1);
                     end
                        
                        
                        
                    if at ~= 0
                        config1(N-nn1+1) = dd1;
                        reclist(nn1-1,dd1,config1,at-1);
                    end
                   
                   
                end
             end
               
          
        
      
             
            reclist(N,1,zeros(2*N,1),at0);
         
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
            nperms=size(permlist,1);
            
            v = permlist;
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
        
        function d = MakeDeltaMatrix(obj,level)
        
        N=obj.nParticles;
        dim = obj.hilbDim;

        sitelist=obj.siteOccupationList(:,1:N);
        atomlist=logical(obj.siteOccupationList(:,(N+1):2*N));
        dlist=zeros(dim,1);
        
        
        for ii = 1:dim
            
            
             v=histc(sitelist(ii,:).*atomlist(ii,:),1:obj.nSites);
             
             
            dlist(ii) = sum(v==level);
        
        end
        
        
        i1 = find(dlist);
        
          
        d = sparse(i1,i1,dlist(i1),dim,dim);   
       

            
        end
        
        function jc = MakeJCMatrix(obj,level)
        
            
            
        a = sparse(kron(eye(obj.nSites),[0,0;1,0]));
        jc=obj.MakeNParticleMatrix(a);
        
        [row col val] = find(jc);
        
        
        jclist = zeros(length(row),3);
        
            
        N=obj.nParticles;
        atoms = N+1:2*N;
        
        
        for ii=1:length(row)
        
            diff = logical(obj.siteOccupationList(row(ii),atoms) - ...
                   obj.siteOccupationList(col(ii),atoms));
               
            
            site = obj.siteOccupationList(row(ii),...
                diff);
               
            %so we know the site at which the transition takes place. Jsut
            %need to work out which level was occupied. This is equal to
            %the number of atoms occu[ied in row(ii) at the site.
            
            if sum(...
                obj.siteOccupationList(row(ii),...
                logical(obj.siteOccupationList(row(ii),atoms)))==site)...
                == level
                jclist(ii,:) = [row(ii),col(ii),val(ii)];
            end
            
            
        end
            
            m = size(jc);
            
        jc = sparse(nonzeros(jclist(:,1)),nonzeros(jclist(:,2)),...
            nonzeros(jclist(:,3)),m(1),m(2));
      
        jc = (jc+jc')/sqrt(level);
        
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
        
            b=[1,0;0,0];
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
                    my=kron(mx,b);
                    mx=kron(my,b);
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
    
     function w = MakeOmegaMatrix(obj)
        
            N = obj.nParticles;
           atoms = obj.siteOccupationList(:,N+1:2*N);
 
           plist=zeros(obj.hilbDim,1);
           for ii =1:obj.hilbDim
           plist(ii)=nnz(obj.siteOccupationList(atoms(ii,:)));
           end
           [ind val] = find(plist);
           w=sparse(ind,ind,val,obj.hilbDim,obj.hilbDim);
     end
         
        
     
     function p = MakePhotonicProjector(obj)
        
           N = obj.nParticles;
           atoms = obj.siteOccupationList(:,N+1:2*N);
           photons=obj.siteOccupationList(logical(atoms));
           
           plist=zeros(obj.hilbDim);
           for ii =1:obj.hilbDim
           plist(ii)=nnz(obj.siteOccupationList(atoms(ii,:)));
           end
           
     end
     
     
     
     function d = MakeDipoleMatrix(obj)
       
        
         [ mx my ] = obj.MakeHardAdjacencyMatrices;
         [ mtx mty ] = obj.MakeWrapAdjacencyMatrices;
         aMat = mx + my + mtx + mty;
         aMat=aMat + aMat' - 2*diag(diag(aMat));
         
        N=obj.nParticles;
        nSites=obj.nSites;
        dim = obj.hilbDim;
       
        
        sitelist=obj.siteOccupationList(:,1:N);
        atomlist=logical(obj.siteOccupationList(:,(N+1):2*N));
        dlist=zeros(dim,1);
        
        for ii = 1:dim
            for jj = 1:N-1
                for kk = jj+1:N
           
                dlist(ii) =dlist(ii) + aMat(sitelist(ii,jj),sitelist(ii,kk))*...
                    atomlist(ii,jj)*atomlist(ii,kk);
                    
                end
            end
        end
        
        
         
        i1 = find(dlist);
        
          
        d = sparse(i1,i1,dlist(i1),dim,dim);   
       
        
         
     end
     
    end %methods
    
   

    
end %classdef

