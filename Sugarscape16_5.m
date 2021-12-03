%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis 930413-6030                                   %%
%% Exercise 16.5 A simple sugarscape                            %%
clear all;

N = 50;
A = 400;
agents = zeros(3,A);

%% Initialize vision, metabolic rate and initial sugar endowment for each 
%% agent
for i = 1:A
    agents(1,i) = unidrnd(6);
    agents(2,i) = unidrnd(4);
    agents(3,i) = unidrnd(21) + 4;
end
vision = agents(1,:);
metabolism = agents(2,:);
sugar = agents(3,:);

%% Initialize sugarscape
 sugarscape = zeros(N,N,2);
% init = 36;
% for i = 1:N/2
%     sugarscape(init,i,1)
%     
% end

%% Randomly place agents in sugarscape
for i = 1:A
    x = randi(N);
    y = randi(N);
    while sugarscape(x,y,2) ~= 0
        x = randi(N);
        y = randi(N);
    end
    sugarscape(x,y,2) = i;
end

%% Simulate movement for each round
rounds = 1e5;
for i = 1:rounds
%% Simulate movement for every agent
for nu = 1:A
agent = randi(A);
[r,c] = find(sugarscape(:,:,2) == agent);
step = vision(agent);
maxSugar = 0;
% check upwards
if (r - step) == 0
    if sugarscape(N,c,1) > maxSugar
        maxSugar = sugarscape(N,c,1);
    end
elseif (r - step) < 1
    excess = abs(r - step);
    for row = (N - excess):N
        if sugarscape(row,c,1) > maxSugar
            maxSugar = sugarscape(row,c,1);
        end
    end
    for row = 1:(r - 1)
        if sugarscape(row,c,1) > maxSugar
            maxSugar = sugarscape(row,c,1);
        end
    end
else
    excess = r - step;
    for row = excess:(r - 1)
        if sugarscape(row,c,1) > maxSugar
            maxSugar = sugarscape(row,c,1);
        end
    end
end
% check downwards
if (r + step) > N
    excess = r + step - N;
    for row = 1:excess
        if sugarscape(row,c,1) > maxSugar
            maxSugar = sugarscape(row,c,1);
        end
    end
    for row = (r + 1):N
        if sugarscape(row,c,1) > maxSugar
            maxSugar = sugarscape(row,c,1);
        end
    end
else
    excess = r + step;
    for row = (r + 1):excess
        if sugarscape(row,c,1) > maxSugar
            maxSugar = sugarscape(row,c,1);
        end
    end
end
% check left
if (c - step) == 0
    if sugarscape(r,column,1) > maxSugar
        maxSugar = sugarscape(r,column,1);
    end
elseif (c - step) < 1
    excess = abs(c - step);
    for column = 1:(c - 1)
        if sugarscape(r,column,1) > maxSugar
            maxSugar = sugarscape(r,column,1);
        end
    end
    for column = (N - excess + 1):N
        if sugarscape(r,column,1) > maxSugar
            maxSugar = sugarscape(r,column,1);
        end
    end
else
    excess = column - step;
    for column = excess:(column - 1)
        if sugarscape(r,column,1) > maxSugar
            maxSugar = sugarscape(r,column,1);
        end
    end
end
% check right     REMAINS TO BE WRITTEN^^
if (c - step) == 0
    if sugarscape(r,column,1) > maxSugar
        maxSugar = sugarscape(r,column,1);
    end
elseif (c - step) < 1
    excess = abs(c - step);
    for column = 1:(c - 1)
        if sugarscape(r,column,1) > maxSugar
            maxSugar = sugarscape(r,column,1);
        end
    end
    for column = (N - excess + 1):N
        if sugarscape(r,column,1) > maxSugar
            maxSugar = sugarscape(r,column,1);
        end
    end
else
    excess = column - step;
    for column = excess:(column - 1)
        if sugarscape(r,column,1) > maxSugar
            maxSugar = sugarscape(r,column,1);
        end
    end
end
        
    





end
end














