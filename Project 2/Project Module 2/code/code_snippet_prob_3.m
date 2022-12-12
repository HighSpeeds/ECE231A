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
  for n_=2:n(length(n))
    % disp(n_)
    new_channel_ps=zeros(1,2*length(channel_ps));
    % disp(size(new_channel_ps))
    for j=1:length(channel_ps)
      new_channel_ps(2*j-1)=channel_ps(j)^2;
      new_channel_ps(2*j)=2*channel_ps(j)-channel_ps(j)^2;
    end
    channel_capacities=ones(1,length(new_channel_ps))-new_channel_ps;
    if any(n==n_)
      for j=1:length(epsilon)
        % size(channel_capacities)
        % sum(channel_capacities<epsilon(j))
        % channel_capacities<epsilon(j)
        output0(n_-4,j)=sum(channel_capacities<=epsilon(j))/length(channel_capacities);
      end
      for j=1:length(epsilon)
        output1(n_-4,j)=sum(channel_capacities>=1-epsilon(j))/length(channel_capacities);
      end
    end
    channel_ps=new_channel_ps;
    % disp("-------------------------------")
  end

  %I move the plotting inside of the loop because its better that way, so then
  %you can see the plots for each p   
  f1 = figure;
  for j=1:length(epsilon)
      plot(n,output1(:,j),'-o','DisplayName',['epsilon = ' num2str(epsilon(j))]);
      hold on;
  end
  grid on;
  legend;
  title(['p=' num2str(p(i))]  );
  xlabel('n');
  ylabel('Fraction of polarized channels with capacity near 1');
  saveas(f1,['p=' num2str(p(i)) ' near 1.png']);

  f2 = figure;
  for j=1:length(epsilon)
      plot(n,output0(:,j),'-o','DisplayName',['epsilon = ' num2str(epsilon(j))]);
      hold on;
  end
  grid on;
  legend;
  xlabel('n');
  title(['p=' num2str(p(i))]  );
  ylabel('Fraction of polarized channels with capacity near 0');
  saveas(f2,['p=' num2str(p(i)) ' near 0.png']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure;
for i=1:length(epsilon)
    plot(n,output1(:,i),'-o','DisplayName',['epsilon = ' num2str(epsilon(i))]);
    hold on;
end
grid on;
legend;
xlabel('n');
ylabel('Fraction of polarized channels with capacity near 1');

f2 = figure;
for i=1:length(epsilon)
    plot(n,output0(:,i),'-o','DisplayName',['epsilon = ' num2str(epsilon(i))]);
    hold on;
end
grid on;
legend;
xlabel('n');
ylabel('Fraction of polarized channels with capacity near 0');

