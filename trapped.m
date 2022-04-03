%% Simulation of Complex Systems                          %%
%% Kakkos Ioannis                                         %%
%% Exercise 6.1 Diffusion regimes for a Brownian particle %%
%% Trapped Brownian Particle                              %%

particles = 10^4;
totalSimTime = 10^1;                    % 10 seconds total sim time
dt = 10^-5;                             % time increment
timesteps = totalSimTime/dt;            % timestep
kB = 1.38064*10^-23;                    % Boltzmann constant
%gamma = 1.88*10^(-8);                   % gas constant
gamma = 1.3;
temperature = 300;                      % standard temperature
trajectory = zeros(10^2,2);
termB = dt/gamma; 
trajectory(1,1) = 0;
trajectory(1,2) = 0;
kx = 10^(-9);
ky = 0.25*10^(-9);
termC = sqrt(2*kB*temperature*dt/gamma);


for i = 2:timesteps
    trajectory(i,1) = trajectory(i-1,1) - termB*kx*trajectory(i-1,1)...
        + termC*randn;
    trajectory(i,2) = trajectory(i-1,2) - termB*ky*trajectory(i-1,2)...
        + termC*randn;
end

% Plotting Inertial Brownian particle
X = trajectory(:,1);
Y = trajectory(:,2);
figure;
hold on
title('Trapped Brownian particle');
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
