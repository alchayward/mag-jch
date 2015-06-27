classdef jch_onsite
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        offset = 100;
        max_n = 20;
        onsiteStrength = struct('b',1,'D',0,'mu',-.78,'k',1,'w',0);
        b
        s
        h0
        dim
        f_psi
    end
    
    methods
       
    function obj = jch_onsite(para)

        if nargin > 0
                obj.onsiteStrength=para.onsiteStrength;
        end
        obj.dim=2*obj.max_n;
        obj.b=obj.MakeBMatrix;
        obj.s=obj.MakeSMatrix;
    end
        
    function h = jchH(obj,h0,psi)
        a = psi'*obj.b;
        h = h0  - ( a + a'); 
        % Note that this is not including the psi*psi terms, since it
        % only shifts the energy, and the interation method does not
        % require it. 
    end
       
    function fPsi = Make_fPsi(obj)
        %First, obtain a polynomial fit of the ground state 
        %order parameter.
        obj.h0 = obj.MakeH0();
        fpsi = @(x)obj.compute_order_paramter( obj.jchH( obj.h0, x) );
        fPsi = @(x)arrayfun(fpsi,x);
    end

    function fPsi = Make_fPsi_fitted(obj,varargin)
      %First, obtain a polynomial fit of the ground state 
      %order parameter.
      n_points = 100; %number of points.
      n_range = [-10,1];
      
      if nargin > 1
        n_points = varargin{1};
      end
      if nargin > 2
         n_range = varargin{2};
      end

      
      
      
      fPsi = obj.Make_fPsi();
      x = [0,logspace(n_range(1),n_range(2),n_points)];
      y = fPsi(x).';
      
      %interp1 is really slow. Handroll it.
      dx = (n_range(2)-n_range(1))/n_points;
      
      
      %fpsi = @(xi)exp(1i*angle(xi))*interp1(x,y,abs(xi));
      fpsi = @(xi)exp(1i*angle(xi))*my_interp1(x,y,abs(xi),dx,n_range);
      fPsi = @(x)arrayfun(fpsi,x);
   end

   function b = MakeBMatrix(obj)
        li = sqrt(1:(obj.max_n-1));
        b=diag(li,1);
        b=kron(b,eye(2));
        b=full(sparse(b));
    end
        
     
   function s =  MakeSMatrix(obj)
        s=kron(eye(obj.max_n),[0,1;0,0]);
        s=full(sparse(s));
   end
       
          
   function h = MakeH0(obj)

        h=obj.onsiteStrength.b.*obj.s'*obj.b;
        h=h+h';
        h=h+obj.onsiteStrength.D.*obj.s'*obj.s;
        h=h-(obj.b'*obj.b + obj.s'*obj.s).*obj.onsiteStrength.mu;

        h=full(sparse(h));
    end
         

   function [a d]= compute_order_paramter(obj,h)
       [v e] = eigs(h,1,'SA');
       gs=v(:,1);
       d=e(1,1);
       %gs=obj.GetGroundstate(h);
       a = gs'*obj.b*gs;
       d=d+obj.offset;
    end

    function V = GetGroundstate(obj,h)
        [V,~]=eigs(h-...
         sparse(1:obj.dim,1:obj.dim,1)*obj.offset,1);   
    end
    end   
    
end

