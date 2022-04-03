%% Simulation of Complex Systems                          %%
%% Kakkos Ioannis                                         %%
%% Exercise 6.4 Continuous time random walk               %%

particles = 10^4;
totalSimTime = 10^-5;                   % 10 microseconds total sim time
dt = 10^-9;                             % time increment
timesteps = 500;                         % timestep
kB = 1.38064*10^-23;                    % Boltzmann constant
gamma = 1.3;                              % gas constant
temperature = 300;                      % standard temperature
%T = zeros(1,totalSimTime);
% trajectory = zeros(timesteps,2);
trajectory = zeros(timesteps,2);
tau = randn*(10^-8);                     % τ = m/gamma a few nanoseconds
termA = (2 + dt/tau)/(1 + dt/tau);
termB = 1/(1 + dt/tau);
termC = (sqrt(2*kB*temperature*gamma)*(dt^(3/2)))/((1 + dt/tau)*gamma*tau);
trajectory(1,1) = 0;
trajectory(2,1) = 0;
trajectory(1,2) = 0;
trajectory(2,2) = 0;

for i = 3:timesteps
    trajectory(i,1) = termA*trajectory(i-1,1) - termB*trajectory(i-2,1)...
        + termC*randn;
    trajectory(i,2) = termA*trajectory(i-1,2) - termB*trajectory(i-2,2)...
        + termC*randn;
end

% Generate waiting times
waitingTimes = zeros(1,timesteps);
r = rand;
alpha = 1;              % for α = 0.5 and α = 1
dt(1) = r^(-1/alpha);
waitingTimes(1) = dt(1);
for i = 2:timesteps
    r = rand;
    dt(i) = r^(-1/alpha);
end
for i = 1:timesteps
    increment = 0;
    for w = 1:i
        increment = increment + dt(w);
    end
    waitingTimes(i) = increment;
end

%% Regularize
[t,P2] = regularize(trajectory,waitingTimes); 

%% Plots
% figure;
% hold on;
% pbaspect([5 2 1]);
% X = waitingTimes;
% Y1 = P1;
% Y2 = P2;
% title('CTRW trajectories for different alpha values');
% xlabel('t[steps]');
% ylabel('x');
% plot(X,Y1,'b','LineWidth',0.8);
% plot(X,Y2,'r','LineWidth',0.8);
% hold off;
% legend('α = 0.5','α = 1');

%% Regularization function
function [Ttime,x] = regularize(irreg,T)

Lirreg = length(irreg);
increments = 2*Lirreg;
time = zeros(1,increments);
tau_next = 0;
tau_prev = 0;
u_next = 0;
u_prev = 0;
for i = 1:Lirreg
    s = 2*i;
    time(s) = T(i);
end
for i = 2:2:increments
    time(i-1) = time(i) - 0.5;
end
for i = 1:increments
    if i <= length(T)
    for j = 1:i
        if T(j) <= time(i)
            if T(j) > tau_prev
                tau_prev = T(j);
            end
            u_prev = randn;
        else
            if T(j) > tau_next
                tau_next = T(j);
            end
            u_next = randn;
        end
    end
    x(i) = u_prev + (u_next - u_prev)*(time(i) - u_prev)/(tau_next - tau_prev);
    end
end

for i = 1:2:length(time) - 1
    Ttime(i) = (time(i) + time(i + 1))/2;
end
Ttime = nonzeros(Ttime);
Ttime = Ttime';
end

