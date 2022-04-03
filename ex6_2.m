%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis                                               %%
%% Exercise 6.2 Regularizing an irregularly sampled trajectory  %%
clear all;clc;

timesteps = 100;
steps = 0;
i = 1;
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
end

s = 1;
I = zeros(1,timesteps);
for j = 1:length(irregular)
    if irregular(j) ~= 0 && s <= timesteps
        I(s) = irregular(j);
        s = s + 1;
    end
end

g = length(irregular);
T = zeros(1,g);
for i = 1:g
    T(i) = i;
end

%% Call function for regularization
[time,regular] = regularize(irregular,T);

%% Call function for normalization
norm = normalise(regular);

%% Plots
figure;
hold on;
grid on;
pbaspect([5 2 1]);
title('Irregularly sampled trajectory');
xlabel('t[steps]');
ylabel('x');
X1 = T;
X2 = time;
Y1 = irregular;
Y2 = regular;
plot(X1,Y1,'o');
plot(X2,Y2);
hold off;
legend('Irregular trajectory','Regularized trajectory');

figure;
hold on;
grid on;
pbaspect([5 2 1]);
title('Normalization');
xlabel('t[steps]');
ylabel('x');
X1 = T;
X2 = time;
Y1 = regular;
Y2 = norm;
plot(X1,Y1,'b');
plot(X2,Y2,'r');
hold off;
legend('Generic trajectory','Normalized trajectory');

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
    x(i) = u_prev + (u_next - u_prev)*(time(i) - tau_prev)/(tau_next - tau_prev);
    end
end

for i = 1:2:length(time) - 1
    Ttime(i) = (time(i) + time(i + 1))/2;
end
Ttime = nonzeros(Ttime);
Ttime = Ttime';
end

function norm_dx = normalise(irreg)
for i = 2:length(irreg)
    dx(i) = irreg(i) - irreg(i - 1);
end
norm_dx = (dx - min(dx))/(max(dx) - min(dx));

end




