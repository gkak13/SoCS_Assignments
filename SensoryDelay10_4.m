%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis                                               %%
%% Exercise 10.4 Multiple robots with sensory delay             %%
 clear all; clc;

n = 50;
L = 10;
tau = 2; % sec
r_node = 1;
i_node = 5;
v_node = 1;
v_inf = 1e-3;

% randomly placing the robots
robots = zeros(10);
c = 0;
while c <= 50
    x = randi(10);
    y = randi(10);
    if (robots(x,y) == 0 && x ~= y)
        robots(x,y) = 1;
        c = c + 1;
    end
end

timesteps = 10^3;

for t = 1:timesteps
s = 0;
for i = 1:L
    for j = 1:L
        if i ~= j
            if robots(i,j) == 1
                s = s + 1;
                pos1(s) = i;
                pos2(s) = j;
                if rem(s,2) == 0
                    dif1(s-1) = (pos1(s) - pos1(s-1))^2;
                    dif2(s-1) = (pos2(s) - pos2(s-1))^2;
                end
            end
        end
    end
end

int(t) = i_node*exp(-(sum(dif1) + sum(dif2))/r_node^2);
end




















