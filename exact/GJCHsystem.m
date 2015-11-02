classdef GJCHsystem < Latticesystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

   % energies=struct('kappa',1,'delta',0,'jc',1);
    angle
    atomLevels = 2;
    end
    
    methods
     function obj = GJCHsystem(para)
            
            obj = obj@Latticesystem(para); 
            if nargin > 0
            obj.atomLevels = para.atomLevels;       
            end               
            obj.model = 'GJCH';         
     end
     
	function p = Initilize(p)
         p.hilb_dim = p.MakeHilbDim;
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

        function k = hopping(obj)
           k = obj.hoppingstrength(1)*obj.mat.k;
        end

        function d = MakeHilbDim(obj)
            d=0;
            np=obj.nParticles;
            p=partitions(np);
            sites=prod(obj.lattice_dim);
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
        
        function LIST = MakeSiteOccupationList(obj)
	%Revisiting this function a year later. I have no f'ing idea how it works.
	%Please don't break it.
         np=obj.nParticles;
         n_sites=obj.nSites;
         at0=obj.atomLevels-1;
         global LIST
         LIST = zeros(obj.hilb_dim,2*np);
         global II
         II=1;
             
            function reclist(np1,site1,config1,at1)
		% np1: particles left to place
		% site1: site that we're up to
		% config1: the configuration we're up to
		% at:	the atomic states of the particle we're up to.
                if np1 == 0 %no more paricles to place
                    LIST(II,:) = config1; % add a new configurations
                    II = II + 1; %next configuration
                else
                    for site = site1:n_sites
                        config1(np-np1+1) = site;
                        reclist(np1-1,site,config1,0);
                    end
                    config1(2*np-np1+1) = 1;
                    for site = (site1+1):n_sites
                        config1(np-np1+1) = site;
                        reclist(np1-1,site,config1,at0-1);
                    end
                    if at1 ~= 0
                        config1(np-np1+1) = site1;
                        reclist(np1-1,site1,config1,at1-1);
                    end
                end
             end
          reclist(np,1,zeros(2*np,1),at0);
        end
                     
        function sym = MakeSymMatrix(obj)
            %for the JCH model, the basis structure is a little more
            %complicated. But we follow a similar proceedure.
            %The basis for distinguishable particles will be 
            %BoxAoxBoxA....Np times. which will be 2^Np times as large as
            %the BH hilb_dim.
            
            len=factorial(obj.nParticles);
            %posm maps a list of particle positions onto a basis index.
            
            %store occuaption as (sss...,aaa...);
            posm=(ones(len,1)*((obj.nParticles-1):-1:0));
            posmb=2*(2*obj.nSites).^posm;
            posma=(2*obj.nSites).^posm;
            
            posm=[posmb';posma']';
            list=zeros(obj.hilb_dim*len,3);
            iter=1;
            
            %%% The calls to perms is very expensive, so we do it by hand.
            %%% permlist is the permutations to loop over;
            
            permlistb  = perms(1:obj.nParticles);
            permlista  = perms(1+obj.nParticles:2*obj.nParticles);
            permlist = [permlistb,permlista];
            nperms=size(permlist,1);
            
            v = permlist;
            shiftsites = [-ones(1,obj.nParticles),zeros(1,obj.nParticles)];
            
            for ii=1:obj.hilb_dim
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
            obj.hilb_dim,(2*obj.nSites)^obj.nParticles);
        end                  
        
        function d = MakeDeltaMatrix(obj,level)
        
        N=obj.nParticles;
        dim = obj.hilb_dim;

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
              
        function [tx ty ] = twist(obj,tangle) %Why?
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
                    [ mx,my ] = obj.MakeHardAdjacencyMatrices;
                    mtx = [];
                    mty = [];
                    mx=kron(mx,b);
                    my=kron(my,b);
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
                    n=round(abs(obj.alpha)*prod(obj.lattice_dim));
                    e=HarperMin(prod(obj.lattice_dim));
                    k=e(n,2);
                  end
                  ssite=[0,beta;...
                      beta,delta-...
                      kappa*k];
                  
                  [v,~]=eig(ssite); 
                  fac=v(:,1);
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

        function PsiL = PfaffianState(obj,varargin)
            
            params = struct();
            twist = obj.tangle;
            if nargin > 1
                params = varargin{1};
                if isfield(params,'twist')
                    twist = params.twist;
                end
            end
            facs = obj.atomic_state_factor(1:obj.hilb_dim,params);
            np = obj.nParticles;
            Z = reshape(...
                obj.sitePositions(obj.siteOccupationList(:,1:np)',1)+...
                obj.sitePositions(obj.siteOccupationList(:,1:np)',2)*1i,...
                np,obj.hilb_dim).';
            lattice_dims = obj.lattice_dim;
            
            wf = Wavefunctions();
            PsiL = wf.Subspace('pfaffian',Z,lattice_dims,twist);
            dim_Psi = size(PsiL,2);
            PsiL = PsiL.*(facs*ones(1,dim_Psi));
            h = HubbardLibrary();
            PsiL = h.GramSchmit(PsiL);
        end
        
        function pro = PNProjector(obj,pn)
	%Not sure what this does now either
           p=obj.nParticles;
           pro=zeros(1,obj.hilb_dim);
           if pn ==1
            for ii = 1:obj.hilb_dim
               if length(unique(nonzeros(obj.siteOccupationList(ii,1:p))))==p
                   pro(ii)=1;
               end
            end
             pro=diag(sparse(pro));
           elseif pn==2
               pro=sparse(1:obj.hilb_dim,1:obj.hilb_dim,1)-obj.PNProjector(1);
           elseif pn==3
           for ii = 1:obj.hilb_dim
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
    
     function w = MakeOmegaMatrix(obj) %I wish i knew what this did.
        
            N = obj.nParticles;
           atoms = obj.siteOccupationList(:,N+1:2*N);
 
           plist=zeros(obj.hilb_dim,1);
           for ii =1:obj.hilb_dim
           plist(ii)=nnz(obj.siteOccupationList(atoms(ii,:)));
           end
           [ind,val] = find(plist);
           w=sparse(ind,ind,val,obj.hilb_dim,obj.hilb_dim);
     end     
     
     function p = MakePhotonicProjector(obj)
        
           N = obj.nParticles;
           atoms = obj.siteOccupationList(:,N+1:2*N);
           photons=obj.siteOccupationList(logical(atoms));
           
           plist=zeros(obj.hilb_dim);
           for ii =1:obj.hilb_dim
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
        dim = obj.hilb_dim;
       
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
