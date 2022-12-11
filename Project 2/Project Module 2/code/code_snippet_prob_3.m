%% Code snippets for ECE 231A: Information Theory: Porject Module #3
% Problem 3 code snippets
clear all;
clc;
close all;

p = [0.3,0.6,0.8]; % Different erasure probabilities for BEC
epsilon = [0.1,0.01,0.001,0.0001];
n = 5:25; % log2 block length

% output for capacity near 0
% rows: corresponding to different n and
% columns: corresponding to different epsilon
output0 = zeros(length(n),length(epsilon));
% output for capacity
output1 = zeros(length(n),length(epsilon));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Enter your code here %%%%%%%%
for i=1:length(p)
  channel_ps=[p(i)^2,2*p(i)-p(i)^2];
  channel_fractions=zeros(n(length(n)));
  for n_=2:n(length(n))
    new_channel_ps=zeros(2*length(channel_ps));
    for j=1:length(channel_ps)
      new_channel_ps(2*j-1)=channel_ps(j)^2;
      new_channel_ps(2*j)=2*channel_ps(j)-channel_ps(j)^2;
    end
    channel_capacities=ones(length(new_channel_ps));
    channel_capacities=channel_capacities.-new_channel_ps;
    if any(n==n_)
      for j=1:length(epsilon)
        sum(channel_capacities<epsilon(j))
        channel_capacities<epsilon(j)
        output_0(n_-4,j)=sum(channel_capacities<epsilon(j))/length(channel_capacities)
      end
      for j=1:length(epsilon)
        output_1(n_-4,j)=sum(channel_capacities<1-epsilon(j))/length(channel_capacities)
      end
    end
    channel_ps=new_channel_ps;
  end
  #channel_ps
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
##f1 = figure;
##for i=1:length(epsilon)
##    plot(n,output1(:,i),'-o','DisplayName',['epsilon = ' num2str(epsilon(i))]);
##    hold on;
##end
##grid on;
##legend;
##xlabel('n');
##ylabel('Fraction of polarized channels with capacity near 1');
##
##f2 = figure;
##for i=1:length(epsilon)
##    plot(n,output0(:,i),'-o','DisplayName',['epsilon = ' num2str(epsilon(i))]);
##    hold on;
##end
##grid on;
##legend;
##xlabel('n');
##ylabel('Fraction of polarized channels with capacity near 0');

