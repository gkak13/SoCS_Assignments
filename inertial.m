%% Simulation of Complex Systems                          %%
%% Kakkos Ioannis 930413-6030                             %%
%% Exercise 6.1 Diffusion regimes for a Brownian particle %%
%% Inertial Brownian particle                             %%

particles = 10^4;
totalSimTime = 10^-5;                   % 10 microseconds total sim time
dt = 10^-9;                             % time increment
timesteps = totalSimTime/dt;            % timestep
kB = 1.38064*10^-23;                    % Boltzmann constant
gamma = 1.3;                              % gas constant
temperature = 300;                      % standard temperature
%T = zeros(1,totalSimTime);
% trajectory = zeros(timesteps,2);
trajectory = zeros(10^4,2);
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

% Plotting Inertial Brownian particle
X = trajectory(:,1);
Y = trajectory(:,2);
figure;
hold on
title('Inertial Brownian particle');
xlabel('x[m]');
ylabel('y[m]');
plot(X,Y);
hold off

% MSD function calling
Tmsd = (0.5:0.5:10);

[eMSDx,eMSDy] = eMSD(trajectory);
%[tMSDx,tMSDy] = tMSD(trajectory);

% Plotting ensemble and time-averaged MSD
figure;
hold on;
title('Ensemble MSDx');
xlabel('t[s]');
ylabel('MSD[m^2]');
plot(Tmsd,eMSDx,'o');
%plot(Tmsd,tMSDx,'*');
hold off;
%legend('ensemble','time-averaged');

% figure;
% hold on;
% title('Ensemble MSDy');
% xlabel('t[s]');
% ylabel('MSD[m^2]');
% plot(Tmsd,eMSDy,'o');
% plot(Tmsd,tMSDy,'*');
% hold off;
% legend('ensemble','time-averaged');

function [eMSDx,eMSDy] = eMSD(path)
mother = (500:500:10^4);
tau = 2;
n = length(path);
ensmbl = zeros(n);
eMSDx = zeros(n);
eMSDy = zeros(n);
step = 1;
RHSx = 0;
RHSy = 0;
for i = 2:n
    RHSx = RHSx + (path(i,1) - path(1,1))^2;
    RHSy = RHSy + (path(i,2) - path(1,2))^2;
    if i == mother(step)
        step = step + 1;
        ex(step) = RHSx;
        ey(step) = RHSy;
    end
end
eMSDx = ex(2:end);
eMSDy = ey(2:end);
end