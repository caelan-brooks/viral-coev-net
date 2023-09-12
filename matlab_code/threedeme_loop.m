clear; close all; clc;

%editing code - nice

% mrates = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3];
% mrates = [1e-8, 1e-7, 5e-7, 1e-6, 1e-5, 1e-4, 5e-4, 1e-3];
mrates = [1e-5];
probs = zeros(length(mrates),1);
probvar = zeros(length(mrates),1);
for mit = 1:length(mrates)
mit
R = 50;
Ncomp = zeros(R,1);
obreak =  zeros(R,1);


for iter = 1:R
% Model parameters
alpha = 0; beta = 1.0; gamma = 1; Fmax = 0.05; D = 1e-3; M = 8; Nh = 1e10; r = 7;
T = 5000;
L = 100; 
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
m12 = mrates(mit);
m23 = mrates(mit);
m13 = mrates(mit);

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


Nkeep = zeros((T*0.1),1);
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

    % migration  
%     for pop = 1:NumPop
%         n_other_pop = mod(pop, 2) + 1;  % the other population
%         n(pop,:) = n(pop,:) + migration_rate * n(n_other_pop,:) - migration_rate * n(pop,:);
%     end

    n(1,:) = n(1,:) + m12*n(2,:) + m13*n(3,:) - m13 * n(1,:) - m12 * n(1,:);
    n(2,:) = n(2,:) + m12*n(1,:) + m23*n(3,:) - m12 * n(2,:) - m23 * n(2,:);
    n(3,:) = n(3,:) + m13*n(1,:) + m23*n(2,:) - m13 * n(3,:) - m23 * n(3,:);


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


    if ii>=(T*0.9)
        Nkeep(ii-(T*0.9)+1) = sum(n(1,:)) + sum(n(2,:)) + sum(n(3,:));
    end 

end
toc;


Ncomp(iter) = mean(Nkeep);
if mean(Nkeep)>1e3
    obreak(iter) = 1;
else 
    obreak(iter) = 0;
end 


end 

probs(mit) = mean(obreak); 
p = mean(obreak);
probvar(mit) = p*(1-p)/(T*0.1);

end 


plot(log(mrates), probs, 'O-')
p2 = probs;
probvar2 = p2.*(1-p2)/(T*0.1);

