
%% define constants
lambda1 = 4 ;
lambda2 = 1 ;
a = 1 ;
k = 60 ;


%% solve odes
initial_condition_odes = zeros(1,k) ;
initial_condition_odes(6) = 1000 ;

[t_meanfield_q12,X_meanfield_q12] = ode15s(@(t,X)RHS_meanfield_q10(t,X,lambda1,a,k),[0,1000],initial_condition) ;   
initial_condition = zeros(k,1) ;

%% use ode solution at t = 50 as I.C for stochastic simulation

% number of repeats
M = 100 ; 

% final time simulate to
T_final = 50 ; 

% recording time interval
rec_step = 0.01 ; 

% vector of time points
time_vec = 0:rec_step:T_final ;

% number of steps required
num_rec_steps = T_final / rec_step ; 

% mean number of particles matrix
mean_num_parts = zeros(k, num_rec_steps + 1);

% Stochastic simulation

% run through a for loop for each repeat
for m = 1:M
    % define initial time
    t = 0 ;
    
    % initialise times to help recording
    t_after = t ;
   
    % define vector for number of particles in each state
    num_parts = zeros(k, 1) ;
   
    % define initial condition (all start in state 0)
    num_parts(1) = initial_condition_odes ;
    
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
                a_f(i,j) = (f(i-j,a)/N)*num_parts(i)*num_parts(j) ;
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

%% make the plots
% find indexes of time points 10 and 50
ind_10 = 10 / rec_step + 1 ; 
ind_50 = 50 /rec_step + 1 ;

states_vector = 0:k-1 ;
f1 = figure ;
f2 = figure ;

% plot bar graphs of mean number of particles in each state
% t = 10
figure(f1)
bar(states_vector, mean_num_parts(:,ind_10))

% t = 50
figure(f2)
bar(states_vector, mean_num_parts(:,ind_50))

%% define f(z) = -a*z*H(-z)
function f = f(z,a)
if z < 0
    f = -a*z ;
else
    f = 0 ;
end
end

%% define mean-field equations from q10
function dX = RHS_meanfield_q10(t,X,lambda1,lambda2,a,k,N)

dX(1) = a.*X(1).*khat(X,k,N) - lambda1.*X(1) + lambda2.*X(2);

for i = 2:k-1
    dX(i) = a.*X(i).*(khat(X,k,N) - (i-1)) + lambda1.*X(i-1) - lambda1.*X(i) ;
end

dX(k) = a.*X(k).*(khat(X,k,N) - (k-1)) + lambda1.*X(k-1) ;

dX = dX' ;
end

function khat = khat(X,k)
   cumsum = 0 ;
   for i = 1:k
       cumsum = cumsum + (i-1).*X(i) ;
   end
   khat = cumsum / sum(X) ; 
end