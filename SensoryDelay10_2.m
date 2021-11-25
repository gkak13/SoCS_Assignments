%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 10.2 Robot in a Gaussian light intensity            %%
% clear all; clc;

tau = 1; 
v0 = 1;                         % max speed [m/s]
v_Inf = 1e-2;                   % degenerate speed [m/s]
timesteps = 10^3;
x = zeros(timesteps,1);         
y = zeros(timesteps,1);   
x1 = zeros(timesteps,1);         
y1 = zeros(timesteps,1); 
x2 = zeros(timesteps,1);         
y2 = zeros(timesteps,1); 
phi = 0;
phi1 = 0;
phi2 = 0;
i_node = 1;
rho_node = 1;
delta = 0;                       % no sensory delay
delta1 = 5*tau;                  % experimenting with different δ values
delta2 = - 5*tau;                % experimenting with different δ values
int = zeros(timesteps,1);
int1 = zeros(timesteps,1);
int2 = zeros(timesteps,1);
index = (timesteps/20:timesteps/20:timesteps);
%index = (5:5:timesteps);
i = 1;
r = zeros(20,1);
r1 = zeros(20,1);
r2 = zeros(20,1);
dist = zeros(timesteps/200,1);
dist1 = zeros(timesteps/200,1);
dist2 = zeros(timesteps/200,1);

for t = 5:timesteps
    v = v_Inf + (v0 - v_Inf)*exp(-int(t)*(t - delta));
    v1 = v_Inf + (v0 - v_Inf)*exp(-int(t)*(t - delta1));
    v2 = v_Inf + (v0 - v_Inf)*exp(-int(t)*(t - delta2));
    white = randn;
    white1 = randn;
    white2 = randn;
    phi = phi + sqrt(2/tau)*white; 
    phi1 = phi1 + sqrt(2/tau)*white1; 
    phi2 = phi2 + sqrt(2/tau)*white2; 
    x(t) = x(t) + v*cos(phi);
    x1(t) = x1(t) + v1*cos(phi1);
    x2(t) = x2(t) + v2*cos(phi2);
    y(t) = y(t) + v*sin(phi); 
    y1(t) = y1(t) + v1*sin(phi1); 
    y2(t) = y2(t) + v2*sin(phi2); 
    int(t) = i_node*exp(-(x(t)+y(t))^2/rho_node^2);
    int1(t) = i_node*exp(-(x1(t)+y1(t))^2/rho_node^2);
    int2(t) = i_node*exp(-(x2(t)+y2(t))^2/rho_node^2);
    dist(t) = sqrt(x(t)^2 + y(t)^2);
    dist1(t) = sqrt(x1(t)^2 + y1(t)^2);
    dist2(t) = sqrt(x2(t)^2 + y2(t)^2);
    if t == index(i)
        r(i) = mean(dist);
        r1(i) = mean(dist1);
        r2(i) = mean(dist2);
        i = i + 1;
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
text(0.8 , 0.8 , txt);
hold off;

figure;
hold on;
title('Mean distance from source point','Interpreter','Latex');
xlabel('t','Interpreter','Latex');
ylabel('$$\langle r(t) \rangle$$','Interpreter','Latex');
plot(r);
plot(r1,'r');
plot(r2,'green');
hold off;
legend('δ = 0','δ = 5τ','δ = -5τ');

