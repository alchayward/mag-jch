function f = thetaln(z,tau,varargin)
% Maps the genrealized Jacobi elliptical theta function to theta3. Returns
% the log of theta[a,b](z|tau)
iter_lim = 7; % Controls error on the function. Will work fine for most cases, but is inefficient.

if nargin == 3 %take in a standard elliptic theta (1,2,3,4)
    switch varargin{1}
        case 1
            %theta_1
            a = 1/2;
            b = 1/2;
        case 2
            %theta_2
            a = 1/2;
            b = 0;
        case 3
            %theta_3
            a = 0;
            b = 0;
        case 4
            %theta_4
            a = 0;
            b = 1/2;
    end
elseif nargin == 4
    a = varargin{1};
    b = varargin{2};
end
    
% Translate to the first lattice square (0,0),(1,1), which makes convergence much better(predicatable). 
a0 = mod(a,1);
z1=z+b+a0*tau;
f_mod = 1i*pi*(tau*a0^2+2*a0*(z+b)); 

%%%
%%% I have put the code from theta3_ln in here, because matlab has terrible 
%%% function calling overhead, and I need to call this a lot. 
%%%

%remove periodic elements so that z lies in the (0,0) to (1,tau) reigion.
%This makes gaurentees about convergence better. Assuming tau = 1 means
%iter_lim = 7 is very well converged. Not true for 'extreme' values of tau

n = floor(imag(z1)./imag(tau));
z2 = z1 - n*tau;
m = floor(real(z2));
z3 = z2 - m;

iter = -iter_lim:iter_lim;
f_sum = sum(exp(1i*pi*tau*ones(length(z),1)*(iter.^2)+2i*pi*z3*iter),2);

f = -1i*pi*(tau*(n.^2)+2*n.*z3)+log(f_sum)+f_mod;
end
