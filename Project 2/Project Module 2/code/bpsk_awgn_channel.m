function [y,Ly] = bpsk_awgn_channel(x,sigma)
%%%% This function is a BPSK + AWGN channel
%%%% x is the input to the channel - Polar encoded data
%%%% sigma is the Gaussian noise std
%%%% The function returns the output from the channel y
%%%% and the corresponding log likelihood ratio Ly
N = length(x);
s = 2*x-1;
y  = s + sigma * randn(1,N); 
Ly = -2*y/sigma^2;
return