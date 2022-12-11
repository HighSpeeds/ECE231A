%% Code snippets for ECE 231A: Information Theory: Porject Module #3
% Problem 4 code snippets
clear all;
close all;
clc;
warning off;

% scenario
    SNR_db = 0:2:10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write your code here that computes the achievable %%
    %%% rates for BPSK+AWGN channel at different SNR(in dB)%
    %%% Plot the figure: achievable rate vs SNR (in dB) %%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    max_trials = 1E4;
    seed1 = 0;
    seed2 = 0;
    n = 10; % log2 block length
    R = 1/2; % rate
    N = 2^n;
    K = round(N*R);
    % Construct code by using Binary Erasure Channel heuristic
    z = 1-R; % erasure rate (z is a real number)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Enter your code here from Problem 3 here %%%%%%%%
    %%%%%% that computes the erasure probabilities  %%%%%%%%
    %%%%%% of N = 2^n polarized channels            %%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Sorting the N = 2^n erasure probabilities %%%%%
    [vals, inds] = sort(z,'ascend'); % z is a vector of length N=2^n here
    freeInds = sort(inds(1:K),'ascend'); % Choosing first K small erasure probability channels
    frozInds = sort(inds(K+1:N),'ascend');
    FI = zeros(1,N);
    FV = zeros(1,N);
    FI(frozInds) = 1;
            
    %% initialize statistics
    FER = zeros(1,length(SNR_db)); % block error rate
    BER = zeros(1,length(SNR_db)); % bit error rate
        %% begin test    
            for snrpoint = 1:length(SNR_db)
                fprintf('\nN: %d, Rate: %4.3f, SNR_db = %3.2f\n', N, R, SNR_db(snrpoint));
                %%%% Write your code here %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                rand('seed',seed1);
                randn('seed',seed2);
                % statistics 
                fer =0;
                ber = 0;
                u = zeros(1,N);
                for trial = 1:max_trials
                    % Put information bits in the K coordinates of u
                    u(freeInds) = round(rand(1,K));
                    % Encodes the message bits u to input data to the
                    % channel x
                    x = polar_encode(u);
                    % Modulate and transmit the data through the channel
                    [y,Ly] = bpsk_awgn_channel(x,sigma);
                    % Decode the polar code from the observed y
                    uh = polar_decode(Ly,FI,FV);
                    
                    %%%%% Enter your code here %%%%%%
                    %%%%% Compute the BER and FER %%%
                    
                    
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if mod(trial,1000) == 0
                        fprintf('N: %d, K: %d, SNR = %4.2f, FER: %5.3e, BER: %5.3e, Trials %d\n',N, K,SNR_db(snrpoint),FER(snrpoint), BER(snrpoint),trial);
                    end
                end   
            end
                
    % plot results
   
    %semilogy(SNR_db,BER,'-db')
    %semilogy(SNR_db,FER,'-db')
    
