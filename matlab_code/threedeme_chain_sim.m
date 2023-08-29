clear; close all; clc;

% Model parameters
alpha = 0; beta = 1.2; gamma = 1; Fmax = 0.05; D = 1e-3; M = 8; Nh = 1e10; r = 7;
T = 5000;
L = 200; 
x = -L:L;
Nx = length(x);
nc = 1;

% Number of populations
NumPop = 3;

n = zeros(NumPop, Nx);
h = zeros(NumPop, Nx);
n_hist = zeros(T, NumPop, Nx);
h_hist = zeros(T, NumPop, Nx);
c_hist = zeros(T, NumPop, Nx);
N_hist = zeros(T, NumPop);

% migration rate
m12 = 1e-7;
m23 = 1e-7;

for pop = 1
    n(pop,:) = exp(-x.^2/5); n(pop,:) = n(pop,:) / sum(n(pop,:)) * Nh / 100;
%     h(pop,:) = heaviside(x+40).*heaviside(-x);
%     h(pop,:) = h(pop,:) / sum(h(pop,:));
    h(pop,:) = 0;
end

for pop = 2
    n(pop,:) = 0;
%     h(pop,:) = heaviside(x+40).*heaviside(-x);
%     h(pop,:) = h(pop,:) / sum(h(pop,:));
    h(pop,:) = 0;
end

for pop = 3
    n(pop,:) = 0;
%     h(pop,:) = heaviside(x+40).*heaviside(-x);
%     h(pop,:) = h(pop,:) / sum(h(pop,:));
    h(pop,:) = 0;
end

tic;

 

for ii = 1:T
    
    c = zeros(NumPop, Nx);
    for pop = 1:NumPop
        for i = 1:Nx
            for k = 1:Nx
                diff_index = min(abs(i - k),Nx - abs(i - k));
                c(pop,i) = c(pop,i) + exp(-diff_index / r) * h(pop,k);
            end
        end

        S = (1 - c(pop,:)).^M;
        F = Fmax * (beta * S - alpha - gamma);
        
        % mutation and growth
        n(pop,:) = n(pop,:) + D * (circshift(n(pop,:),1) + circshift(n(pop,:),-1) - 2 * n(pop,:));
        n(pop,:) = (1 + F) .* n(pop,:);
    end


    n(1,:) = n(1,:) + m12*n(2,:) - m12 * n(1,:);
    n(2,:) = n(2,:) + m12*n(1,:) + m23*n(3,:) - m12 * n(2,:) - m23 * n(2,:);
    n(3,:) = n(3,:) + m23*n(2,:) - m23 * n(3,:);


    for pop = 1:NumPop
        N = sum(n(pop,:));
        h(pop,:) = h(pop,:) + 1/M/Nh * (n(pop,:)- N.*h(pop,:));

        n(pop,:) = poissrnd(n(pop,:));        


    
        n_hist(ii,pop,:) = n(pop,:);
        h_hist(ii,pop,:) = h(pop,:);
        c_hist(ii,pop,:) = c(pop,:);
        N_hist(ii,pop) = N;
        
        F_hist(ii,:) = F;

      
        
    end



end
toc;

figure;
for j=1:10:T
    for pop = 1:NumPop
        subplot(3,1,pop);
        S = (1 - squeeze(c_hist(j,pop,:))) .^ M; F = (beta * S - alpha - gamma);

        
        plot(x, squeeze(n_hist(j,pop,:)), 'r','LineWidth',2);
        yyaxis left
%         hold on
        plot(x, squeeze(h_hist(j,pop,:)), 'b','LineWidth',2);  
        yyaxis right



%         plot(x, squeeze(F_hist(j,:)), 'b','LineWidth',2);  
%         yyaxis right
%         plot(x, Nh * squeeze(h_hist(j,pop,:)), 'b','LineWidth',2);
%         plot(x, F, 'k','LineWidth',2);
%         plot(x, squeeze(c_hist(j,pop,:)), 'k','LineWidth',2); hold off;
%         hold off
        
%         title(sprintf('Population %d', pop));
%         xlabel('x');
%         ylabel('Density');
       
        
        
%         ylim([0 max(n_hist(:))])
    end

    drawnow;
    pause(0.05);
end

