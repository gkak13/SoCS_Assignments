%% Simulation of Complex Systems                          %%
%% Kakkos Ioannis                                         %%
%% Exercise 6.1 Diffusion regimes for a Brownian particle %%
%% Overdamped Brownian particle                           %%

particles = 10^4;
totalSimTime = 10^1;                    % 10 seconds total sim time
trajectory = zeros(particles,2);
dt = 10^-3;                             % time increment
timesteps = totalSimTime/dt;            % timestep
kB = 1.38064*10^-23;                    % Boltzmann constant
%gamma = 1.3;                            % gas constant
gamma = 1.88*10^(-8); 
temperature = 300;                      % standard temperature
term = sqrt(2*kB*temperature*dt/gamma);
T = zeros(1,totalSimTime);
eMSDx = zeros(timesteps,1);
eMSDy = zeros(timesteps,1);

% create timestep array
for j = 1:timesteps - 1
    T(j) = j*dt;
end
T(10^4) = 10^4*dt;

previousX = 0;
previousY = 0;
for i = 1:particles
    % x-direction
    r = randn;
    nextX = previousX + term * r;
    previousX = nextX;
    trajectory(i,1) = nextX;
   
    % y-direction
    q = randn;
    nextY = previousY + term * q;
    previousY = nextY;
    trajectory(i,2) = nextY;
end

% create time array for eMSD-tMSD
Tmsd = (0.5:0.5:10);

[eMSDx,eMSDy] = eMSD(trajectory);
% [tMSDx,tMSDy] = tMSD(trajectory);

% Plotting overdamped Brownian particle
X = trajectory(:,1);
Y = trajectory(:,2);
figure;
hold on
title('Overdamped Brownian particle');
xlabel('x[m]');
ylabel('y[m]');
plot(X,Y);
hold off

% Plotting ensemble and time-averaged MSD
figure;
hold on;
title('Ensemble MSDx');
xlabel('t[s]');
ylabel('MSD[m^2]');
plot(Tmsd,eMSDx,'o');
% plot(Tmsd,tMSDx,'*');
hold off;
% legend('ensemble','time-averaged');

figure;
hold on;
title('Ensemble MSDy');
xlabel('t[s]');
ylabel('MSD[m^2]');
plot(Tmsd,eMSDy,'o');
% plot(Tmsd,tMSDy,'*');
hold off;
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

% function [tMSDx,tMSDy] = tMSD(path)
% mother = (500:500:10^4);
% tau = 2;
% n = length(path);
% Times = zeros(n);
% tMSDx = zeros(n);
% tMSDy = zeros(n);
% step = 1;
% RHSx = 0;
% RHSy = 0;
% for i = 2:n
%     RHSx = RHSx + (path(i,1) - path(1,1))^2;
%     RHSy = RHSy + (path(i,2) - path(1,2))^2;
%     if i == mother(step)
%         step = step + 1;
%         tx(step) = RHSx;
%         ty(step) = RHSy;
%     end
% end
% tMSDx = tx(2:end);
% tMSDy = ty(2:end);
% end
