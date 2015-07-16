function f = theta(a,b,z,tau)

iterMax = 200;
inds = -iterMax:iterMax;
theta_fun = @(n)exp(1i*pi*tau*(a+n)^2+2i*pi*(n+a)*(z+b));
f = sum(arrayfun(theta_fun,inds));
end