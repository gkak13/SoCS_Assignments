%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis                                               %%
%% Exercise 10.3 Robot in a circular well                       %%
% clear all; clc;

tau = 1; v0 = 1; v_Inf = 1e-1; radius = 2;
realizations = 1e5;
timesteps = 1e3;
delta = 5*tau;
x = zeros(timesteps,1); y = zeros(timesteps,1);    
phi = 0; i_node = 1; rho_node = 1;                       
int = zeros(timesteps,1);

for i = 1:realizations
for t = 6:timesteps
    int(t) = i_node*exp(-(x(t)+y(t))^2/rho_node^2);
    if delta > 0
        v = v_Inf + (v0 - v_Inf)*exp(-int(t - delta));
    else
        diff = int(t) - int(t - 3);
        v = v_Inf + (v0 - v_Inf)*exp(-(int(t) - delta*diff*int(t)));
    end
    white = randn;
    phi = phi + sqrt(2/tau)*white;   
    x(t + 1) = x(t) + v*cos(phi);
    y(t + 1) = y(t) + v*sin(phi); 
    dist(t) = x(t + 1)^2 + y(t + 1)^2;
    r(t) = sqrt(x(t)^2 + y(t)^2);
    if dist(t) > radius^2
        x(t + 1) = x(t)*radius/sqrt(x(t)^2 + y(t)^2);
        y(t + 1) = y(t)*radius/sqrt(x(t)^2 + y(t)^2);
    end
end
rf(i) = mean(r);                 % first stage averaging
end

s = 1;
for i = 1:realizations
    if rem(i,10) == 0
        rs(s) = mean(rf(1:i));    % second stage averaging                                               
        s = s + 1;
    end
end

figure;
hold on;
title('Robot motion - Circular well, negative delay','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x,y);
txt = {'Data:',['δ = ',num2str(delta)],['Io = ',num2str(i_node)],...
    ['ρο = ',num2str(rho_node)]};
text(max(x), max(y), txt);
hold off;

figure;
hold on;
title('Radial drift - negative delay','Interpreter','Latex');
xlabel('$$x^{2} + y^{2}$$','Interpreter','Latex');
ylabel('$$\langle r(t) \rangle$$','Interpreter','Latex');
plot(rs);
hold off;
