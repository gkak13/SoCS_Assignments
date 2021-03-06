%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis                                               %%
%% Exercise 10.1 Simulation of a light sensitive robot          %%
% c.
% clear all; clc;

tau = 1; 
c = 1e-5;       % aribtrary choice for constant c
v0 = 5*10^-6;   % max speed [m/s]
v_Inf = 10^-6;  % degenerate speed [m/s]
timesteps = 10^4;
lamda = 5e-6;
x = zeros(timesteps,1);
y = zeros(timesteps,1);
phi = 0;
int = (sin(2*pi*1e-6/lamda))^2;   % arbitrary choice that x(0) = 1e-6

for t = 1:timesteps
    int = (sin(2*pi*(x(t)-c*t)/lamda))^2;
    v = v_Inf + (v0 - v_Inf)*exp(-int);
    white = randn;
    phi = phi + sqrt(2/tau)*white; 
    x(t + 1) = x(t) + v*cos(phi);
    y(t + 1) = y(t) + v*sin(phi); 
end

% plot
figure;
hold on;
title('Robot motion - Time-varying light pattern','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x,y);
txt = {'Data:',['cτ = ',num2str(c)],['Λ = ',num2str(lamda)],...
    ['L = 5e-06']};
text(max(x), max(y), txt);
hold off;









































