%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 10.2 Robot in a Gaussian light intensity            %%
clear all; clc;

realizations = 1e4;
timesteps = 1e2;
tau = 1; 
v0 = 1;                         
v_Inf = 1e-1;                   
% x = zeros(timesteps,1); y = zeros(timesteps,1);   
% x1 = zeros(timesteps,1); y1 = zeros(timesteps,1); 
% x2 = zeros(timesteps,1); y2 = zeros(timesteps,1); 
phi = 0; phi1 = 0; phi2 = 0;
i_node = 1; rho_node = 1;
delta = 0; delta1 = 5*tau; delta2 = -5*tau;                 
int = zeros(timesteps,1); int1 = zeros(timesteps,1);
int2 = zeros(timesteps,1);

for i = 1:realizations
    x = zeros(timesteps,1); y = zeros(timesteps,1);   
    x1 = zeros(timesteps,1); y1 = zeros(timesteps,1); 
    x2 = zeros(timesteps,1); y2 = zeros(timesteps,1); 
for t = 2:timesteps
    int(t) = i_node*exp(-(x(t)+y(t))^2/rho_node^2);
    int1(t) = i_node*exp(-(x1(t)+y1(t))^2/rho_node^2);
    int2(t) = i_node*exp(-(x2(t)+y2(t))^2/rho_node^2);
    if delta > 0
        v = v_Inf + (v0 - v_Inf)*exp(-int(t - delta));
    else
        diff = int(t) - int(t - 1);
        v = v_Inf + (v0 - v_Inf)*exp(-(int(t) - delta*diff*int(t-1)));
    end
    diff = int(t) - int(t - 1);
    v1 = v_Inf + (v0 - v_Inf)*exp(-(int(t) - delta1*...
        diff*int(t-1)));  % case of δ > 0
    v2 = v_Inf + (v0 - v_Inf)*exp(-(int(t) - delta2*...
        diff*int(t-1)));  % case of δ < 0
    white = randn;
    white1 = randn;
    white2 = randn;
    phi = phi + sqrt(2/tau)*white; 
    x(t + 1) = x(t) + v*cos(phi);
    x1(t + 1) = x1(t) + v1*cos(phi);
    x2(t + 1) = x2(t) + v2*cos(phi);
    y(t + 1) = y(t) + v*sin(phi); 
    y1(t + 1) = y1(t) + v1*sin(phi); 
    y2(t + 1) = y2(t) + v2*sin(phi); 
%     int(t) = i_node*exp(-(x(t)+y(t))^2/rho_node^2);
%     int1(t) = i_node*exp(-(x1(t)+y1(t))^2/rho_node^2);
%     int2(t) = i_node*exp(-(x2(t)+y2(t))^2/rho_node^2);
    dist(t) = sqrt(x(t)^2 + y(t)^2);
    dist1(t) = sqrt(x1(t)^2 + y1(t)^2);
    dist2(t) = sqrt(x2(t)^2 + y2(t)^2);
end
r(i) = mean(dist);
r1(i) = mean(dist1);
r2(i) = mean(dist2);
end

% figure;
% hold on;
% title('Robot motion - No sensory delay','Interpreter','Latex');
% xlabel('x[m]','Interpreter','Latex');
% ylabel('y[m]','Interpreter','Latex');
% plot(x,y);
% txt = {'Data:',['δ = ',num2str(delta)],['Io = ',num2str(i_node)],...
%     ['ρο = ',num2str(rho_node)]};
% text(max(x), max(y), txt);
% hold off;
% 
% figure;
% hold on;
% title('Robot motion - Sensory delay','Interpreter','Latex');
% xlabel('x[m]','Interpreter','Latex');
% ylabel('y[m]','Interpreter','Latex');
% plot(x1,y1);
% txt = {'Data:',['δ = ',num2str(delta1)],['Io = ',num2str(i_node)],...
%     ['ρο = ',num2str(rho_node)]};
% text(max(x1), max(y1), txt);
% hold off;
% 
% figure;
% hold on;
% title('Robot motion - Sensory delay','Interpreter','Latex');
% xlabel('x[m]','Interpreter','Latex');
% ylabel('y[m]','Interpreter','Latex');
% plot(x2,y2);
% txt = {'Data:',['δ = ',num2str(delta2)],['Io = ',num2str(i_node)],...
%     ['ρο = ',num2str(rho_node)]};
% text(max(x2), max(y2), txt);
% hold off;

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
