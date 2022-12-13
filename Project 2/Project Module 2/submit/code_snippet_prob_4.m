%% Code snippets for ECE 231A: Information Theory: Porject Module #3
% Problem 4 code snippets
clear all;
close all;
clc;
warning off;

% scenario
    SNR_db = 0:2:10;
    size(SNR_db)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write your code here that computes the achievable %%
    %%% rates for BPSK+AWGN channel at different SNR(in dB)%
    %%% Plot the figure: achievable rate vs SNR (in dB) %%%%
    
    capacities=zeros(1,length(SNR_db));

    for i=1:length(SNR_db)
        sigma = sqrt(10^(-SNR_db(i)/10));
        f = @(t) 1/(sigma*sqrt(2*pi))*exp(-(t-1).^2/(2*sigma^2)).*log2(1+exp(-2*t/(sigma^2)));
        capcity = 1-integral(f,-20,20);
        capacities(i) = capcity;
    end
    
    figure(1)
    plot(SNR_db,capacities,'-db')
    xlabel('SNR (dB)')
    ylabel('Achievable Rate (bits/symbol)')
    title('Achievable Rate vs SNR (dB) for BPSK+AWGN channel')
    grid on;
    saveas(gcf,'Achievable Rate vs SNR (dB) for BPSK+AWGN channel.png')

FERs=zeros(5,length(SNR_db));
BERs=zeros(5,length(SNR_db));
Rs=0.1:0.2:0.9;
max_trials = 2E3;
seed1 = 0;
seed2 = 0;
n = 10; % log2 block length
for i=1:length(Rs)
    for snrpoint = 1:length(SNR_db)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R=Rs(i)*capacities(snrpoint);
        N = 2^n;
        K = round(N*R);
        % Construct code by using Binary Erasure Channel heuristic
        z = 1-R; % erasure rate (z is a real number)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Enter your code here from Problem 3 here %%%%%%%%
        %%%%%% that computes the erasure probabilities  %%%%%%%%
        %%%%%% of N = 2^n polarized channels            %%%%%%%%
        
        channel_ps=[2*z - z^2,z^2];
        for n_=2:n(length(n))
            % disp(n_)
            new_channel_ps=zeros(1,2*length(channel_ps));
            % disp(size(new_channel_ps))
            for j=1:length(channel_ps)
                new_channel_ps(2*j-1)=2*channel_ps(j)-channel_ps(j)^2;
                          new_channel_ps(2*j)=channel_ps(j)^2;
    
            end
            channel_ps=new_channel_ps;
            % disp("-------------------------------")
        end
        z = channel_ps;
        for n_ = 0:N-1
          z(n_+1) = channel_ps(bin2dec(fliplr(dec2bin(n_, n)))+1);
        end
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
        fprintf('\nN: %d, Rate: %4.3f, SNR_db = %3.2f\n', N, R, SNR_db(snrpoint));
        %%%% Write your code here %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sigma = sqrt(10^(-SNR_db(snrpoint)/10));
        
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
            
            BERs(i,snrpoint) = BERs(i,snrpoint)+sum(uh ~= u);
            FERs(i,snrpoint) = FERs(i,snrpoint)+any(uh ~= u);

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mod(trial,1000) == 0
                fprintf('N: %d, K: %d, SNR = %4.2f, FER: %5.3e, BER: %5.3e, Trials %d\n',N, K,SNR_db(snrpoint),FERs(i,snrpoint)/(N*trial), BERs(i,snrpoint)/trial,trial);
            end
        end   
    end
    BERs=BERs/(N*max_trials);
    FERs=FERs/(max_trials);
    % plot results
end
figure(2)
for snrpoint = 1:length(SNR_db)
    plot(Rs*capacities(snrpoint),BERs(:,snrpoint))
    hold on;
end
xlabel('rate')
ylabel('BER')
title('BER vs rate for diffrent SNR (dB) for BPSK+AWGN channel')
legend(string(SNR_db)+" DB")
grid on;
saveas(gcf,'BER vs rate for diffrent SNR (dB) for BPSK+AWGN channel.png')
figure(3)
for snrpoint = 1:length(SNR_db)
    plot(Rs*capacities(snrpoint),FERs(:,snrpoint))
    hold on;
end
hold on;
xlabel('rate')
ylabel('FER')
title('FER vs rate for diffrent SNR (dB) for BPSK+AWGN channel')
legend(string(SNR_db)+" DB")
grid on;
saveas(gcf,'FER vs rate for diffrent SNR (dB) for BPSK+AWGN channel.png')
    
    
    
