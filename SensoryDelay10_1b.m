%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis                                               %%
%% Exercise 10.1 Simulation of a light sensitive robot          %%
% b.
% clear all; clc;

tau = 1; 
v0 = 5*10^-6;   % max speed [m/s]
v_Inf = 10^-6;  % degenerate speed [m/s]
timesteps = 10^4;
lamda = [10^-6:10^-6:5*10^-6];    % plots with respect to various Λ values
x = zeros(timesteps,length(lamda));
y = zeros(timesteps,length(lamda));
phi = 0;
int = (sin(2*pi*1e-6/lamda(4)))^2;

for j = 1:length(lamda)
for t = 1:timesteps
    int = (sin(2*pi*x(t,j)/lamda(j)))^2;
    v = v_Inf + (v0 - v_Inf)*exp(-int);
    white = randn;
    phi = phi + sqrt(2/tau)*white; 
    x(t + 1,j) = x(t,j) + v*cos(phi);
    y(t + 1,j) = y(t,j) + v*sin(phi); 
end
end

x1 = x(:,1); x2 = x(:,2); x3 = x(:,3); x4 = x(:,4); x5 = x(:,5);
y1 = y(:,1); y2 = y(:,2); y3 = y(:,3); y4 = y(:,4); y5 = y(:,5);

% plots
figure;
hold on;
title('Robot motion - Periodic light condition','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x1,y1);
txt = {'Data:',['τ = ',num2str(tau)],['I = (sin(2πx(t)/Λ))^2'],['Λ = 20% of L']};
text(max(x1), max(y1), txt);
hold off;

figure;
hold on;
title('Robot motion - Periodic light condition','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x2,y2);
txt = {'Data:',['τ = ',num2str(tau)],['I = (sin(2πx(t)/Λ))^2'],['Λ = 40% of L']};
text(max(x2), max(y2), txt);
hold off;

figure;
hold on;
title('Robot motion - Periodic light condition','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x3,y3);
txt = {'Data:',['τ = ',num2str(tau)],['I = (sin(2πx(t)/Λ))^2'],['Λ = 60% of L']};
text(max(x3), max(y3), txt);
hold off;

figure;
hold on;
title('Robot motion - Periodic light condition','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x4,y4);
txt = {'Data:',['τ = ',num2str(tau)],['I = (sin(2πx(t)/Λ))^2'],['Λ = 80% of L']};
text(max(x4), max(y4), txt);
hold off;

figure;
hold on;
title('Robot motion - Periodic light condition','Interpreter','Latex');
xlabel('x[m]','Interpreter','Latex');
ylabel('y[m]','Interpreter','Latex');
plot(x5,y5);
txt = {'Data:',['τ = ',num2str(tau)],['I = (sin(2πx(t)/Λ))^2'],['Λ = 100% of L']};
text(max(x5), max(y5), txt);
hold off;

