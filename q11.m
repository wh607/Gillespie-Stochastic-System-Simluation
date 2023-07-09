%% define constants
% total number of particles
N = 100 ;

a = 1 ;
lambda1 = 5 ;

% number of states
k = 25 ;

% vector of states used for plots later
states_vector = 0:k-1 ;

%% solve q10 odes
initial_condition = 4.* ones(1,k) ;
[t_meanfield, X_meanfield] = ode15s(@(t,X)RHS_meanfield_q10(t,X,lambda1,a,k),[0,1000],initial_condition) ;

index_stable = length(X_meanfield) ;

%% steady state in q9
k_0 = 0 ; % first state in which initial condition non-zero is n_0
q9steadystates = zeros(1,25) ;
% create a vector of theoretical distribution of states given in Question 9
for i = 0:k-1
   q9steadystates(i + 1) = scaledpoisson(N, lambda1, a, i, k_0) ;
end

%% make plots
f1 = figure ;
figure(f1)
plot(states_vector, q9steadystates, 'linewidth', 5)
hold on
plot(states_vector, X_meanfield(index_stable , :), 'Color', 'r', 'LineStyle', '--', 'linewidth' , 3)
xlabel('Particle state')
ylabel('Number of particles')
title('Steady states: Theoretical distribution vs ODE solutions')
legend('Theoretical distribution', 'Solution of ODEs')

%% define scaled poisson eq'n
function steadystate = scaledpoisson(N,lambda1, a, k, k_0)
    steadystate = N*exp(-lambda1/a)*(lambda1/a)^(k-k_0)*(1/factorial(k-k_0)) ;
end

% define mean-field equations from q10
function dX = RHS_meanfield_q10(t,X,lambda1,a,k)

    dX(1) = a.*X(1).*khat(X,k) - lambda1.*X(1) ;

    for i = 2:k-1
        dX(i) = a.*X(i).*(khat(X,k) - (i-1)) + lambda1.*X(i-1) - lambda1.*X(i) ;
    end

    dX(k) = a.*X(k).*(khat(X,k) - (k-1)) + lambda1.*X(k-1) ;

    dX = dX' ;
end

% define function that calculates k^
function khat = khat(X,k)
   cumsum = 0 ;
   for i = 1:k
       cumsum = cumsum + (i-1).*X(i) ;
   end
   khat = cumsum / sum(X) ;
end
