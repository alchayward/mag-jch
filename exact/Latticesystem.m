classdef Latticesystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    model='BH' 
    alpha=0;
    alphap=0;
    alphaq=1;
    A;
    tangle = [0,0];
    nParticles=1;
    hilb_dim;
    lattice_dim = [4,4];
    mat
    nSites = 1;
    sitePositions
    topology = 'torus' %Other option is 'Hard', with no wrapping.
    siteOccupationList;
    sym; 
    bosonList;
    laughlinwf;
    onsiteStrength;
    hoppingStrength;
    pot=0;
    maxParticlesPerSite = 1;
    sysPath = './'
    V;
    D;
    
    %Lattice Properties

        %Presently this code only handles the square lattice geometery.
        
        %Single Particle Adjacency Matrtrix
        %it is nessasary to store seperate X and Y adjacency matricies, as well
        %as X and Y tunneling matricies to make changing the twist angles
        %easier. These are stored as 3-tuples (i,j,c-number)

    %Operators
        %Operators, like the Jaynes-Cummings interaction, and hopping matrix,
        %are stored as sparse matricies, and only with their upper right
        %triangle elements.

        %Occupation Lists

        %These lists are Np by Ns matricies, storing the position and number
        %of photonic and atomic occupations at each site.

    
    end
    
    methods
    
        
        %Constructor
     function obj = Latticesystem(para)
            if nargin > 0
                    obj.lattice_dim = para.lattice_dim;
                    obj.nParticles = para.nParticles;
                    obj.alphap=para.alpha(1);
                    obj.alphaq=para.alpha(2);
                    obj.alpha = obj.alphap/obj.alphaq;
                    obj.onsiteStrength=para.onsiteStrength;
                    obj.hoppingStrength=para.hoppingStrength;
                    obj.A=para.A;
                    obj.maxParticlesPerSite = para.maxParticlesPerSite;
                    obj.sysPath = para.sysPath;
                    obj.tangle=para.twist;
            end

            obj.nSites=obj.NumberOfSites;
            obj.sitePositions=obj.MakeSitePositions;
          %  obj.laughlinwf=laughlin(para);
        end
          
        %System Things
 
        function m = MakeNParticleMatrix(obj,a)
            
        an=a;
        id=length(a);
        for ii=1:(obj.nParticles-1)
           idd=id^ii;
           e = speye(id);
            an=kron(a,sparse(1:idd,1:idd,1)) + kron(e,an);
        end

        %reduce to ind basis;
        m=obj.sym*an*obj.sym';
        

        end      
        
        function sitepos = MakeSitePositions(obj)
            
            %Sites are referenced by a linear index. This fucntion recovers
            %the (x,y) coordinates 
            
            sitepos= zeros(obj.nSites,2);
            [sitepos(:,1),sitepos(:,2)] = ... 
            ind2sub(obj.lattice_dim,1:obj.nSites);
            sitepos(:,1)=sitepos(:,1)-1;
            sitepos(:,2)=sitepos(:,2)-1;
        end  
                      
        function n = NumberOfSites(obj)
           n = obj.lattice_dim(1) * obj.lattice_dim(2);
        end

        function [mx,my] = MakeHardAdjacencyMatrices(obj)
            % gives the Adjacency matrix for the 2D lattice.
            %   Sites are labeled from 1 to Lx*Ly starting bottom left,
            %and working horizontally to top right.
            %

            %The topology defines what the lattice looks like. at the moment, only
            %option is for a square lattice with wither torroidal or hard boundry
            %conditions.

            %dimensios of lattice
            mx = zeros((obj.lattice_dim(1)-1)*obj.lattice_dim(2),3);
            my = zeros((obj.lattice_dim(2)-1)*obj.lattice_dim(1),3);
            %Make the X matrix
            %The hopping is from j to i
            pon=obj.pot;
            poff=1-obj.pot;
            
            incr = 1;
            for ii = 1:(obj.nSites)
               for jj = (ii+1):(obj.nSites)
                  if (((obj.sitePositions(jj,1)-obj.sitePositions(ii,1)) == 1)...
                     && (obj.sitePositions(ii,2) == obj.sitePositions(jj,2)))
                        mx(incr,:) = [ii,jj,exp(1i*...
                            obj.A(1)*2*pi*obj.alpha*...
                            obj.sitePositions(ii,2) - ...
                        (pon*obj.sitePositions(ii,1) + poff)*obj.tangle(1))];
                        incr = incr + 1;
                  end %if

                end %for
             end %for

            %Make Y matrix
            
              
            incr = 1;
            for ii = 1:(obj.nSites)
               for jj = (ii+1):(obj.nSites)
                  if (((obj.sitePositions(jj,2)-obj.sitePositions(ii,2)) == 1)...
                     && (obj.sitePositions(ii,1) == obj.sitePositions(jj,1)))
                        
                        my(incr,:) = [ii,jj,exp(1i*...
                            obj.A(2)*2*pi*obj.alpha*...
                            obj.sitePositions(ii,1) -...
                            (pon*obj.sitePositions(ii,2) + poff)*obj.tangle(2))];
                        incr = incr + 1;
                  end %if

                end %for
             end %for
             
             mx=conj(sparse(mx(:,1),mx(:,2),mx(:,3),obj.nSites,obj.nSites));
             my=conj(sparse(my(:,1),my(:,2),my(:,3),obj.nSites,obj.nSites));
%                  mx=mx+mx';
%                  my=my+my';
        end %Make A matrix
        
        function [mx, my] = MakeWrapAdjacencyMatrices(obj)
           
                
            
            
            % gives the Wrap around Adjacency matrix for the 2D lattice.
            %Sites are labeled from 1 to Lx*Ly starting bottom left,
            %and working horizontally to top right.
            %

            %The topology defines what the lattice looks like. at the moment, only
            %option is for a square lattice with wither torroidal or hard boundry
            %conditions.

            
            mx = zeros(obj.lattice_dim(2),3);
            my = zeros(obj.lattice_dim(1),3);
            pon=obj.pot;
            poff=1-obj.pot;
            
            %Make the X matrix
            %The hopping is from j to i
            
            
            incr = 1;
             for ii = 1:(obj.nSites)
               for jj = 1:(obj.nSites)
                  if ((obj.sitePositions(ii,1) == 0) && ... 
                      (obj.sitePositions(jj,1) == (obj.lattice_dim(1)-1))...
                      && (obj.sitePositions(ii,2) == obj.sitePositions(jj,2)))
                      
                  phase = obj.tangle(1)*(pon*obj.sitePositions(ii,1) + poff)+...
                    2*pi*obj.alpha*obj.sitePositions(ii,2)*...
                        (obj.A(1) + 0*obj.lattice_dim(1));

                        mx(incr,:) = [ii,jj,exp(-1i*phase)];
                        incr = incr + 1;
                  end %if

                end %for
             end %for

            %Make Y matrix
            
             
            incr = 1;
             for ii = 1:(obj.nSites)
               for jj = 1:(obj.nSites)
                  if ((obj.sitePositions(ii,2) == 0) && ... 
                      (obj.sitePositions(jj,2) == (obj.lattice_dim(2)-1))...
                      && (obj.sitePositions(ii,1) == obj.sitePositions(jj,1)))
                 
%                   phase = (obj.tangle(2)...
%                             +2*pi*obj.alpha*...
%                         obj.sitePositions(ii,1)*(obj.lattice_dim(2)));      
                  
                    
                 phase = obj.tangle(2)*(pon*obj.sitePositions(ii,2) + poff)+...
                    2*pi*obj.alpha*obj.sitePositions(ii,1)*...
                        (obj.A(2) + obj.lattice_dim(2));   
                    
                        my(incr,:) = [ii,jj,...
                            exp(-1i*phase)];
                        incr = incr + 1;
                  end %if

                end %for
             end %for
                         
        mx=conj(sparse(mx(:,1),mx(:,2),mx(:,3),obj.nSites,obj.nSites));
        my=conj(sparse(my(:,1),my(:,2),my(:,3),obj.nSites,obj.nSites));
%         mx=mx+mx';
%         my=my+my';
        end    
        
        function [tx, ty] = MakeTwistMatricies(obj)
%             
%             obj.alpha=0;
%             tw=[pi/2 pi/2];
%             obj.tangle = tw;
%             obj.tangle(2)=0;
%             [mx my mtx mty] = obj.MakeAdjacencyMatricies;
%             
%             a=mx+my+mtx+mty;
%             
%             tx = obj.MakeNParticleMatrix(a+a');
%             
%             
%             
%             obj.tangle=tw;
%             obj.tangle(1)=0;
%             [mx my mtx mty] = obj.MakeAdjacencyMatricies;
%             
%             a=mx+my+mtx+mty;
%             
%             ty = obj.MakeNParticleMatrix(a+a');
%         
            %now just need to normalize each element. There's maybe a way faster 
            %way to do this.
            
            
            
            obj.alpha=0;
            tw=[0 0];
            obj.tangle = tw;
            obj.tangle(2)=0;
            [~, ~, mtx, mty] = obj.MakeAdjacencyMatricies; %#ok<MCNPN>
            tx = logical(obj.MakeNParticleMatrix(mtx));
            ty = logical(obj.MakeNParticleMatrix(mty));

%             
%             [r c s] = find(tx);
%             tx=sparse(r,c,s./abs(s),obj.hilbDim,obj.hilbDim);
%              
%             [r c s] = find(ty);
%             ty=sparse(r,c,s./abs(s),obj.hilbDim,obj.hilbDim);
%             
        end
        
        function k = MakeKappaMatrix(obj)
            %obj.tangle = [0 0];
            [mx, my, mtx, mty] = obj.MakeAdjacencyMatricies; %#ok<MCNPN>
            a=mx+my+mtx+mty;
            k = obj.MakeNParticleMatrix(a+a');
           % m = obj.HoppingStrengths(m);            
        end
        
        %%%                         GenerateFileName
%%%
%%%
%%%
%%%  generates the filename for a system defined by systemParameters
%%%
function fileName = GenerateFileName(p)
   fileName=[p.model 'm',num2str(p.maxParticlesPerSite) ...
        'd' num2str(p.lattice_dim(1)) 'x' num2str(p.lattice_dim(2)) ...
        'n' num2str(p.nParticles) 'p' num2str(p.alphap) ...
        'q' num2str(p.alphaq) '.mat']; 
end

function filePath = GenerateFilePath(p) 
    fileName=p.GenerateFileName;
    filePath=fullfile(p.sysPath,fileName);
end
function Z = coordinates(obj)
Z = reshape(...
       obj.sitePositions(obj.siteOccupationList(:,1:np)',1)+...
       obj.sitePositions(obj.siteOccupationList(:,1:np)',2)*1i,...
       np,obj.hilb_dim).';
end
end %methods    
end %classdef

