%% Define the constants
% total number of particles
N = 100 ; 

lambda1 = 1 ;
lambda2 = 1 ;

% 30 states 0 to 29
k = 30 ; 

% vector of states used for plots later
states_vector = 0:k-1 ;

% number of repeats
M = 1000 ; 

% final time simulate to
T_final = 100 ; 

% recording time interval
rec_step = 0.01 ; 

% vector of time points
time_vec = 0:rec_step:T_final ;

% number of steps required
num_rec_steps = T_final / rec_step ; 

% mean number of particles matrix
mean_num_parts = zeros(k, num_rec_steps + 1);

%% Stochastic simulation
% run through a for loop for each repeat
for m = 1:M
    % define initial time
    t = 0 ;
    
    % initialise times to help recording
    t_after = t ;
   
    % define vector for number of particles in each state
    num_parts = zeros(k, 1) ;
   
    % define initial condition (all start in state 0)
    num_parts(1) = N ;
    
    % update mean number of particles at initial time
    mean_num_parts(:,1) = mean_num_parts(:,1) + num_parts ;
    
    % enter the while loop
    while t < T_final
        
        % calculate the propensity functions
        a_left = num_parts.*lambda2 ;
        a_left(1) = 0; % cannot jump left out of first box
        a_right = num_parts.*lambda1 ;
        a_right(k) = 0 ; % cannot jump right out of last box
        % a_f is a matrix where the i-j element represents a particle
        % in state (i - 1) copying the state of a particle in state (j - 1)
        % [-1's due to difference in matlab and question indexing]
        a_f = zeros(k,k) ;
        for i = 1:k
            for j = 1:k
                a_f(i,j) = (abs(i-j)/N)*num_parts(i)*num_parts(j) ;
            end
        end
        
        % calculate the sum of the propensity functions
        a_leftsum = sum(a_left) ;
        a_rightsum = sum(a_right) ;
        a_fsum = sum(sum(a_f)) ;
        a0 = a_leftsum + a_rightsum + a_fsum ;
        
        % determine time for next reaction
        tau = (1/a0)*log(1/rand) ;
        
        % update time
        t = t + tau ;
        
        % generate a random number to determine which reaction occurs
        r1 = a0 * rand ;
        
        % write the current value of num_parts to num_parts_old for
        % recording purposes
        num_parts_old = num_parts ;
        
        % determine which reaction should occur
        
        if r1 < a_leftsum
            % then a left jump has occured
            cumsum = a_left(1) ;
            j = 1;
            % figure out which left jump has happened
            while cumsum < r1
                j = j + 1 ;
                cumsum = cumsum + a_left(j) ;
            end
            % implement left jump from state (j-1) [as matlab indexes from
            % 1, but our states start from n_0]
            num_parts(j) = num_parts(j) - 1 ;
            num_parts(j - 1) = num_parts(j - 1) + 1 ;
            
        elseif r1 < a_leftsum + a_rightsum
            % then a right jump has occured
            cumsum = a_right(1) ;
            j = 1;
            % figure out which right jump has happened
            while cumsum < r1 - a_leftsum
                j = j + 1 ;
                cumsum = cumsum + a_right(j) ;
            end
            % implement right jump from state (j-1) [as matlab indexes from
            % 1, but our states start from n_0]
            num_parts(j) = num_parts(j) - 1 ;
            num_parts(j + 1) = num_parts(j + 1) + 1 ;
        else
            cumsum = a_f(1,1) ;
            i = 1 ;
            j = 1 ;
            while cumsum < r1 - a_leftsum - a_rightsum 
                for j = 1:k
                    for i = 1:k
                        cumsum = cumsum + a_f(i,j) ;
                    end
                end        
            end
            % implement partcle in state (i - 1) copying the state of
            % particle in state (j - 1)
            num_parts(i) = num_parts(i) - 1 ;
            num_parts(j) = num_parts(j) + 1 ;
        end
        
        % calculate the times for recording
        t_before = t_after ;
        t_after = t ;
        
        % calculate the indices of the time step before and the 
        % current time step in terms of recording
        ind_before = ceil((t_before + eps) / rec_step) ;
        ind_after = min(floor(t_after / rec_step),num_rec_steps) ;
        
        % find out how many time-steps to write to
        steps_to_write = ind_after - ind_before + 1 ;
        
        
        if steps_to_write > 0 && t + tau < inf
            mean_num_parts(:, ind_before + 1:ind_after + 1) = mean_num_parts(:, ind_before + 1:ind_after + 1) + (num_parts_old * ones(1, steps_to_write)) ;
        end
            
        
    end


end

% find mean number of particles in each state at each time point
mean_num_parts = mean_num_parts / M ;


%% solve mean-field ODEs
initial_condition = zeros(k, 1) ; 
initial_condition(1) =  N ;
[t_meanfield, X_meanfield] = ode15s(@(t,X)RHS_meanfield_model(t,X,lambda1,lambda2,k),[0,10,30,T_final],initial_condition) ;


%% make the plots
% find indexes of time points 10, 30 and 100
ind_10 = 10 / rec_step + 1 ; 
ind_30 = 30 / rec_step + 1 ;
ind_100 = 100 /rec_step + 1 ;

% create figures
f1 = figure ;
f2 = figure ;
f3 = figure ;


% plot bar graphs of mean number of particles in each state
% t = 10
figure(f1)
bar(states_vector, mean_num_parts(:, ind_10))
hold on
plot(states_vector, X_meanfield(2,:), 'r', 'linewidth', 5)
xlabel('Particle state')
ylabel('Number of particles')
title('Mean-field ODE solutions vs stochastic simulation at t = 10')
hold off

% t = 30
figure(2)
bar(states_vector,mean_num_parts(:, ind_30))
hold on
plot(states_vector, X_meanfield(3,:), 'r', 'linewidth', 5)
xlabel('Particle state')
ylabel('Number of particles')
title('Mean-field ODE solutions vs stochastic simulation at t = 30')
hold off

% t = 100
figure(3)
bar(states_vector,mean_num_parts(:, ind_100))
hold on
plot(states_vector, X_meanfield(4,:), 'r', 'linewidth', 5)
xlabel('Particle state')
ylabel('Number of particles')
title('Mean-field ODE solutions vs stochastic simulation at t = 100')
hold off


%% define mean-field ODE RHS
% f is even so sum over f(k-i) - f(i-k) vanishes and do not include here
function dX = RHS_meanfield_model(t,X,lambda1,lambda2, k)
dX(1) = -lambda1.*X(1) + lambda2.*X(2) ;

for i = 2:k-1
   dX(i) =  lambda1.*X(i-1) - lambda1.*X(i) + lambda2.*X(i+1) - lambda2.*X(i) ;
end

dX(k) = lambda1.*X(k-1) - lambda2.*X(k) ;

dX = dX' ;
end
