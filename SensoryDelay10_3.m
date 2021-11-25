%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 10.3 Robot in a circular well                       %%
% clear all; clc;


tau = 1; 
v0 = 1;                                 % max speed [m/s]
v_Inf = 1e-1;                           % degenerate speed [m/s]
radius = 2;
timesteps = 8e5;
x = zeros(timesteps,1);         
y = zeros(timesteps,1);    
phi = 0;
i_node = 1;
rho_node = 1;
delta = tau;                       
int = zeros(timesteps,1);
index = (5:timesteps);
i = 1;
% r = zeros(20,1);
% dist = zeros(timesteps/200,1);
s = 1;

for t = 5:timesteps
    if delta > 0
        v = v_Inf + (v0 - v_Inf)*exp(-int(t - delta));
    else
        if ((int(t) > int(t - 1)) && (int(t) > int(t - 3)))
            diff = int(t) - int(t - 3);
            v = v_Inf + (v0 - v_Inf)*exp(-(int(t) - delta*diff*int(t)));
        else
            diff = int(t - 3) - int(t);
            v = v_Inf + (v0 - v_Inf)*exp(-(int(t) - delta*diff*int(t)));
        end
    end

    white = randn;
    phi = phi + sqrt(2/tau)*white;   
    x(t + 1) = x(t) + v*cos(phi);
    y(t + 1) = y(t) + v*sin(phi); 
    int(t + 1) = i_node*exp(-(x(t + 1)+y(t + 1))^2/rho_node^2);
    dist(t) = x(t + 1)^2 + y(t + 1)^2;
    r(t) = sqrt(x(t)^2 + y(t)^2);
    if rem(t,50) == 0     
        rt(s) = mean(r);
        s = s + 1;
    end
    if dist(t) > radius^2
        x(t + 1) = x(t)*radius/sqrt(x(t)^2 + y(t)^2);
        y(t + 1) = y(t)*radius/sqrt(x(t)^2 + y(t)^2);
    end
end

figure;
hold on;
title('Robot motion - Circular well, positive delay','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x,y);
txt = {'Data:',['δ = ',num2str(delta)],['Io = ',num2str(i_node)],...
    ['ρο = ',num2str(rho_node)]};
text(max(x), max(y), txt);
hold off;

figure;
hold on;
title('Radial drift','Interpreter','Latex');
xlabel('$$x^{2} + y^{2}$$','Interpreter','Latex');
ylabel('$$\langle r(t) \rangle$$','Interpreter','Latex');
plot(rt);
hold off;
