%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 10.3 Robot in a circular well                       %%
% clear all; clc;


tau = 0.5; 
v0 = 1;                                 % max speed [m/s]
v_Inf = 1e-1;                           % degenerate speed [m/s]
radius = 2;
timesteps = 10^3;
x = zeros(timesteps,1);         
y = zeros(timesteps,1);    
phi = 0;
i_node = 1;
rho_node = 1;
delta = -0.5*tau;                       % no sensory delay
int = zeros(timesteps,1);
index = (5:5:timesteps);
i = 1;
r = zeros(20,1);
dist = zeros(timesteps/200,1);

for t = 2:timesteps
    v = v_Inf + (v0 - v_Inf)*exp(-int(t)*(t - delta));
    white = randn;
    phi = phi + sqrt(2/tau)*white; 
    x(t) = x(t) + v*cos(phi);
    y(t) = y(t) + v*sin(phi); 
    int(t) = i_node*exp(-(x(t)+y(t))^2/rho_node^2);
    dist(t) = x(t)^2 + y(t)^2;
    r(t) = sqrt(x(t)^2 + y(t)^2);
    if dist(t) > radius^2
        x(t) = x(t)*radius/sqrt(x(t)^2 + y(t)^2);
        y(t) = y(t)*radius/sqrt(x(t)^2 + y(t)^2);
    end
end

figure;
hold on;
title('Robot motion - No sensory delay','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x,y);
txt = {'Data:',['δ = ',num2str(delta)],['Io = ',num2str(i_node)],...
    ['ρο = ',num2str(rho_node)]};
text(0.8, 0.8 , txt);
hold off;

figure;
hold on;
title('Radial drift','Interpreter','Latex');
xlabel('$$x^{2} + y^{2}$$','Interpreter','Latex');
ylabel('$$\langle r(t) \rangle$$','Interpreter','Latex');
plot(dist,r);
hold off;
