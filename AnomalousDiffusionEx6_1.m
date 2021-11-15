%% Simulation of Complex Systems                          %%
%% Kakkos Ioannis 930413-6030                             %%
%% Exercise 6.1 Diffusion regimes for a Brownian particle %%
%% Overdamped Brownian particle                           %%

realizations = 10^5;
trajectory = zeros(realizations,2);
dt = 10^-4;             % timestep
kB = 1.38064*10^-23;    % Boltzmann constant
gamma = 1.3;            % gas constant
temperature = 285;      % standard temperature
term = sqrt(2*kB*temperature*dt/gamma);
previousX = 0;
previousY = 0;

for i = 1:realizations
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

ensemble = eMSD(trajectory);
% time = tMSD(trajectory);



function ensmbl = eMSD(path)
RHS = 0;
tau = 1;
for i = 1:length(path)
    RHS = (path(i + tau,1) - path(tau,1))^2;
    ensmbl(i) = RHS; 
end
ensmbl = RHS/length(path);
end

% function t = tMSD(path)
% 
% end


