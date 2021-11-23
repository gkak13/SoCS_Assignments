%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 10.1 Simulation of a light sensitive robot          %%

clear all;clc;

x_t = v*cos(phi_t);
y_t = v*sin(phi_t);
phi_t = sqrt(2/tau)*white;
v = v_Inf + (v0 - v_Inf)*exp(-int);
timesteps = 10^3;


% absence of light
int = 0; % light intensity
for t = 1:timesteps
    x(t) = x(t) + v*cos(phi(tau,white));
    y(t) = y(t) + v*sin(phi(tau,white));
    
    
end

function orientation = phi(tau,w)
orientation = sqrt(2/tau)*w;
end

    









