function list = HarperMin( p,q )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


cl=makemf([p,1],0);

if nargin == 2
    cl.alpha=(q)/p;
    A=cl.MakeAdjacencyMatrix;
    list=min(eig(full(A)));
else
    list=zeros(p-1,2);
    for ii=1:p-1
        cl.alpha=(ii)/p;
        A=cl.MakeAdjacencyMatrix;
        list(ii,:)=[ii/p,min(eig(full(A)))];
    end
end
end

