function f = theta3ln(z,tau)
iter_lim = 7;

%remove periodic elements so that z lies in the (0,0) to (1,tau) reigion.
%This makes gaurentees about convergence better. Assuming tau = 1 means
%iter_lim = 7 is very well converged.
n = floor(imag(z)./imag(tau));
z1 = z - n*tau;
m = floor(real(z1));
z2 = z1 - m;

iter = -iter_lim:iter_lim;
f_sum = sum(exp(1i*pi*tau*ones(length(z),1)*(iter.^2)+2i*pi*z2*iter),2);

f = -1i*pi*(tau*(n.^2)+2*n.*z2)+log(f_sum);
end