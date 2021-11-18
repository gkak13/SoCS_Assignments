%% Simulation of Complex Systems                          %%
%% Kakkos Ioannis 930413-6030                             %%
%% Exercise 6.3 Scaled Brownian Motion                    %%

%clear all;clc;

timesteps = 100;
steps = 0;
i = 1;
alpha = 2;
while steps <= timesteps
    i = i + 1;
    r = randn;
    q = randn;
    if r > 0
        irregular(i) = q;
        steps = steps + 1;
    else
        irregular(i) = nan;
        steps = steps + 1;
    end
    c = randn;
    p = randn;
    if c > 0
        irregular2(i) = p;
        steps2 = steps + 1;
    else
        irregular2(i) = nan;
        steps2 = steps2 + 1;
    end
end

s = 1;
I = zeros(1,timesteps);
for j = 1:length(irregular)
    if irregular(j) ~= 0 && s <= timesteps
        I(s) = irregular(j);
        s = s + 1;
    end
end
s = 1;
I = zeros(1,timesteps);
for j = 1:length(irregular2)
    if irregular2(j) ~= 0 && s <= timesteps
        I(s) = irregular2(j);
        s = s + 1;
    end
end

g = length(irregular);
T = zeros(1,g);
for i = 1:g
    T(i) = i^(1/alpha);
end

%% Call function for regularization
[time,regular] = regularize(irregular,T);
[time,regularX3] = regularize(irregular2,T);

% figure;
% hold on;
% title('2D Plane for different alpha values');
% xlabel('x');
% ylabel('y');
% Y = regular;
% X1 = regularX1;
% X2 = regularX2;
% X3 = regularX3;
% plot(X1,Y,'g','LineWidth',0.8);
% plot(X2,Y,'b','LineWidth',0.8);
% plot(X3,Y,'r','LineWidth',0.8);
% hold off;
% legend('α = 0,5','α = 1','α = 2');

% Run 3 times changing alpha (0.5,1,2) and regular variable name in order
% to store the results for cummulative plotting

% figure;
% hold on;
% grid on;
% pbaspect([5 2 1]);
% title('Scaled Brownian Motion');
% xlabel('t[steps]');
% ylabel('x');
% X = time;
% Y1 = regular;
% Y2 = regular1;
% Y3 = regular2;
% plot(X,Y1,'g','LineWidth',1);
% plot(X,Y2,'b','LineWidth',1);
% plot(X,Y3,'r','LineWidth',1);
% hold off;
% legend('α = 0.5','α = 1','α = 2');

% 2D Plots


%% Regularization function
function [Ttime,x] = regularize(irreg,T)

Lirreg = length(irreg);
increments = 2*Lirreg;
time = zeros(1,increments);
tau_next = 0;
tau_prev = 0;
u_next = 0;
u_prev = 0;
for i = 1:Lirreg
    s = 2*i;
    time(s) = T(i);
end
for i = 2:2:increments
    time(i-1) = time(i) - 0.5;
end
for i = 1:increments
    if i <= length(T)
    for j = 1:i
        if T(j) <= time(i)
            if T(j) > tau_prev
                tau_prev = T(j);
            end
            u_prev = randn;
        else
            if T(j) > tau_next
                tau_next = T(j);
            end
            u_next = randn;
        end
    end
    x(i) = u_prev + (u_next - u_prev)*(time(i) - u_prev)/(tau_next - tau_prev);
    end
end

for i = 1:2:length(time) - 1
    Ttime(i) = (time(i) + time(i + 1))/2;
end
Ttime = nonzeros(Ttime);
Ttime = Ttime';
end
