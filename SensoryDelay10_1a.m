%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 10.1 Simulation of a light sensitive robot          %%
% a.
% clear all; clc;

int = 0;        % light intensity - absence of light
tau = 2; 
v0 = 5*10^-6;   % max speed [m/s]
v_Inf = 10^-6;  % degenerate speed [m/s]
v = v_Inf + (v0 - v_Inf)*exp(-int);
timesteps = 10^5;
x = zeros(1,timesteps);
y = zeros(1,timesteps);
phi = 0;

for t = 1:timesteps
    white = randn;
    phi = phi + sqrt(2/tau)*white; 
    x(t+1) = x(t) + v*cos(phi);
    y(t+1) = y(t) + v*sin(phi); 
end

% plot
figure;
hold on;
title('Robot motion - Absence of light','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x,y);
txt = {'Data:',['τ = ',num2str(tau)],['I = ',num2str(int)]};
text(max(x) , max(y) , txt);
hold off;











