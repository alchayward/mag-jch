function f = theta(a,b,z,tau)

   iterLim=20;

   a0 = mod(a,1);
   z1=z-b;
   
   n = floor(imag(z1)/imag(tau));
   m = floor(real(z1)-n*real(tau));
   z0 = z1 - m - n*tau;
   factor = 1i*pi*(2*a0*m-2*n*(z0+b)-tau*n^2);

    f=0;
    
    for ii=-iterLim:iterLim
        f=f+exp(1i*pi*tau*(ii+a0)^2)*exp(2i*pi*(ii+a)*z0);
    end
        f=f*exp(factor);
%    

%    f=0;
%    for ii=-iterlimit:iterlimit
%       m = exp(1i*pi*(t*(ii+a)^2 +2*(ii+a)*(z+b)));
% 	  f=f+m;
% 								   
% 	 
%     end
end