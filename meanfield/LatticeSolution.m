function f = LatticeSolution(alpha,dim,tau,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ml = MagneticLattice();
l = 1/sqrt(alpha*2*pi);
Lx = dim(1);
Ns = dim(1)*dim(2)*alpha;

f=exp(-imag((z/l))^2/2)/sqrt(Lx/l*Ns*sqrt(pi));
f=f*ml.theta(-0*Ns/Ns,0,z*Ns/Lx,tau*Ns);
if isnan(f)
    f=0;
end

end

