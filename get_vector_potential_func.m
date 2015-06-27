function [ A_func ] = get_vector_potential_func( gauge, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if strcmp(gauge,'landau')
    A_func = @(X)2*pi*alpha*[0;X(1)];
elseif strcmp(gauge,'symmetric')
    A_func = @(X)2*pi*alpha*[-X(2);X(1)]/2;
end
end

