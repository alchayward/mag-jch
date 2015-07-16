classdef BHsystem < Latticesystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties   
    %energies=struct('kappa',1,'delta',0,'jc',1);
    
    end
    
    methods
    
        
        %Constructor
        function obj = BHsystem(para)
            obj = obj@Latticesystem(para); 
            
            if nargin == 0
                   
            end
            
            
             obj.model = 'BH';
             
          
             obj.hilbDim = obj.MakeHilbDim;
        end
          
        function p = Initilize(p)
         
         

         p.siteOccupationList=p.MakeSiteOccupationList;
         
         
         p.sym=p.MakeSymMatrix;

         p.mat.k=p.MakeKappaMatrix;
         
         nBodyInteractions = 1:p.maxParticlesPerSite;
         p.mat.u = {};
         for ii = nBodyInteractions
            p.mat.u{ii} = p.MakeBOnsite(ii);
         end
         
         
        end
        
        
        %System Things
        function   h = MakeHamiltonian(obj,varargin)
                
            nvarargin = length(varargin);
                if nvarargin == 0
                    k=obj.hoppingStrength;
                    u = 0;
                elseif nvarargin == 2
                    k = varargin{1};
                    u = varargin{2};
                else
                    disp(['ERROR - BHystem:MakeHamiltonian:',...
                        'wrong number of arguments']);
                end
                
                h =  - k*obj.mat.k;
                
                oslen = min([length(u),obj.maxParticlesPerSite]);
                
                for ii = 1:oslen
                    h = h + u(ii)*obj.mat.u{ii};
                end
        end
       
  
        function k = hopping(obj)
           k = obj.hoppingstrength*obj.mat.k;
        end
        
        function [tx ty ] = twist(obj,tangle)
            obj.tangle=tangle;
            [tx ty] = obj.MakeTwistMatricies;
            
        end
              
        function d = MakeHilbDim(obj)
            
            
            
            d=0;
            np=obj.nParticles;
            maxN =obj.maxParticlesPerSite;
            p=partitions(np);
            
            sites=prod(obj.latticeDim);
            
            if maxN >= np
                d = nchoosek(obj.nParticles+obj.nSites -1,obj.nParticles);
            else
                
                
                
                
                for ii = 1:size(p,1)
                    pp=p(ii,:);
                    flag = ~nnz(pp(maxN+1:np));
                    if flag
                    fac=factorial(sum(pp))/prod(factorial(pp));
                    else
                        fac = 0;
                    end
                    
                    d=d+fac*nchoosek(sites,sum(pp));
                    
                end
                
            end
            
            
            
            
            
            
            
        end
        
        
        
        
        function sym = MakeSymMatrix(obj)
  
            len=factorial(obj.nParticles);
            
            posm=ones(len,1)*((obj.nParticles-1):-1:0);

            posm=obj.nSites.^posm;
            
            
            
            list=zeros(obj.hilbDim*len,3);
            iter=1;
            for ii=1:obj.hilbDim
                v=perms(obj.siteOccupationList(ii,:));
                li=unique((sum(posm.*(v-1),2)+1)');
                fac=nnz(li);
                nlist=[ones(1,fac)*ii;li;ones(1,fac)/sqrt(fac)]';
                list(iter:iter+fac-1,:)=nlist;
                iter=iter+fac;
            end
            list=list(1:iter-1,:);
            sym=sparse(list(:,1),list(:,2),list(:,3),...
            obj.hilbDim,obj.nSites^obj.nParticles);
         end
                   
        function list = MakeSiteOccupationList(obj,varargin)

         N=obj.nParticles;
         
          maxN =obj.maxParticlesPerSite;
         if length(varargin) == 1
             maxN = varargin{1};
         end
         
         
         d=obj.nSites;
        
             %initilize the list to zero.
             list = zeros(obj.hilbDim,N);
             ii=1;
             function reclist(nn,dd,config,limit)
                if nn == 0
                    list(ii,:) = config;
                    ii = ii + 1;
                else
                    if limit ~= 0
                        config(N-nn+1) = dd;
                        reclist(nn-1,dd,config,limit-1)
                        
                    end
                    
                    for jj = dd+1:d                     
                        config(N-nn+1) = jj;
                        reclist(nn-1,jj,config,maxN-1)
                    end
                end
             end
             
            reclist(N,1,zeros(N,1),maxN);
        
        end
        
        function m = MakeBOnsite(obj,varargin)
        
            
            
            
            nbodies = 2;
            if length(varargin) == 1
                nbodies = varargin{1};
            end
            
            
        m=zeros(obj.hilbDim,1);
        
        for ii=1:obj.hilbDim
             y=sort(obj.siteOccupationList(ii,:));
             v=histc(obj.siteOccupationList(ii,:),y);
            
            u=1;
            for jj = 0:(nbodies-1)
               u =(v-jj).*u;
            end
               u=u/factorial(nbodies);
             m(ii)=sum(u);
            
        end
        
        [row, ~, val] = find(m);
        m=sparse(row,row,val,obj.hilbDim,obj.hilbDim);
        end

        function [mx my mtx mty] = MakeAdjacencyMatricies(obj)
        
             switch obj.topology
                case 'torus'
                    [ mx my ] = obj.MakeHardAdjacencyMatrices;
                    [ mtx mty ] = obj.MakeWrapAdjacencyMatrices;
                case 'hard'
                    [ mx my ] = obj.MakeHardAdjacencyMatrices;
                    mtx = 0;
                    mty = 0;
             end
                    
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
            
     function li = SingleParticlePositions(obj,ind)
        
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
     
     
     function e = E0est(obj)
        
         [p,q] = rat(obj.alpha);
         e=obj.hoppingstrength*obj.nParticles*HarperMin(q,p);
     end
    
     
    function psiL = WavefunctionState(obj,varargin)
       
        
            params = struct();
            twist = obj.tangle;
            wf_type = 'laughlin';
            if nargin > 1
                params = varargin{1};
                if isfield(params,'twist')
                    twist = params.twist;
                end
                if isfield(params,'wf_type')
                    wf_type = params.wf_type;
                end
            end

            Z = reshape(...
                obj.sitePositions(obj.siteOccupationList',1)+...
                obj.sitePositions(obj.siteOccupationList',2)*1i,...
                np,obj.hilbDim).';
            
            lattice_dims = obj.latticeDim;
            
            psiL = zeros(obj.hilbDim,2);
            wf = Wavefunctions(wf_type);
            
            psiL(:,1) = facs.*wf(Z,lattice_dims,1,twist);
            psiL(:,2) = facs.*wf(Z,lattice_dims,2,twist);

            h = HubbardLibrary();
            psiL = h.GramSchmit(psiL);
       
    end
     
    end %methods
    
    
end %classdef

