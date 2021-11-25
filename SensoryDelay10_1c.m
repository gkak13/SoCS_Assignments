%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 10.1 Simulation of a light sensitive robot          %%
% c.
% clear all; clc;

tau = 2; 
c = 1e-6;          % aribtrary choice for constant c
v0 = 5*10^-6;   % max speed [m/s]
v_Inf = 10^-6;  % degenerate speed [m/s]
timesteps = 10^3;
lamda = 1e-6;
x = zeros(timesteps,1);
y = zeros(timesteps,1);
phi = 0;
int = (sin(2*pi*1e-6/lamda))^2;   % arbitrary choice that x(0) = 1e-6

for t = 1:timesteps
    v = v_Inf + (v0 - v_Inf)*exp(-int);
    white = randn;
    phi = phi + sqrt(2/tau)*white; 
    x(t) = x(t) + v*cos(phi);
    y(t) = y(t) + v*sin(phi); 
    int = (sin(2*pi*(x(t)-c*t)/lamda))^2;
end

% plot
figure;
hold on;
title('Robot motion - Time-varying light pattern','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x,y);
txt = {'Data:',['τ = ',num2str(tau)],['I = (sin(2π(x-ct)/Λ))^2']};
text(3.9*10^-6 , 3.9*10^-6 , txt);
hold off;









































