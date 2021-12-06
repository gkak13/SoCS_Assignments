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
x01 = 36; y01 = 13; x02 = 13; y02 = 36;
r1 = 5; r2 = 4; r3 = 4; r4 = 4;
for i = 1:N
for j = 1:N
    dif1 = sqrt((i - x01)^2 + (j - y01)^2);
    dif2 = sqrt((i - x02)^2 + (j - y02)^2);
    if dif1 <= r1
        sugarscape(i,j) = 4;
    elseif (dif1 > r1 && dif1 <= (r1 + r2))
        sugarscape(i,j) = 3;
    elseif (dif1 > (r1 + r2) && dif1 <= (r1 + r2 + r3))
        sugarscape(i,j) = 2;
    elseif (dif1 > (r1 + r2 + r3) && dif1 <= (r1 + r2 + r3 + r4))
        sugarscape(i,j) = 1;
    end
    if dif2 <= r1
        sugarscape(i,j) = 4;
    elseif (dif2 > r1 && dif2 <= (r1 + r2))
        sugarscape(i,j) = 3;
    elseif (dif2 > (r1 + r2) && dif2 <= (r1 + r2 + r3))
        sugarscape(i,j) = 2;
    elseif (dif2 > (r1 + r2 + r3) && dif2 <= (r1 + r2 + r3 + r4))
        sugarscape(i,j) = 1;
    end
end
end

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
rounds = 5e2;
chosen = zeros(1,A);
for i = 1:rounds
%% Simulate movement for every agent
for nu = 1:A
agent = randi(A);
chosen(nu) = agent;
while (find(agent == chosen(1:nu-1)) ~= 0 |...
        find(agent == chosen(nu+1:end)) ~= 0)
    agent = randi(A);
end
chosen(nu) = agent;
[r,c] = find(sugarscape(:,:,2) == agent);
step = vision(agent);
maxSugar = zeros(1,2);
m = 1;
% check upwards
if (r - step) == 0
    if sugarscape(N,c,1) > maxSugar(:)
        maxSugar(m) = sugarscape(N,c,1);
    end
    m = 1;
elseif (r - step) < 1
    excess = abs(r - step);
    for row = (N - excess):N
        if sugarscape(row,c,1) > maxSugar(:)
            maxSugar(m) = sugarscape(row,c,1);
            m = m + 1;
        end
    end
    for row = 1:(r - 1)
        if sugarscape(row,c,1) > maxSugar(:)
            maxSugar(m) = sugarscape(row,c,1);
        end
    end
    m = 1;
else
    excess = r - step;
    for row = excess:(r - 1)
        if sugarscape(row,c,1) > maxSugar(:)
            maxSugar(m) = sugarscape(row,c,1);
        end
    end
    m = 1;
end
% check downwards
if (r + step) > N
    excess = r + step - N;
    for row = 1:excess
        if sugarscape(row,c,1) > maxSugar(:)
            maxSugar(m) = sugarscape(row,c,1);
            m = m + 1;
        end
    end
    for row = (r + 1):N
        if sugarscape(row,c,1) > maxSugar(:)
            maxSugar(m) = sugarscape(row,c,1);
        end
    end
    m = 1;
else
    excess = r + step;
    for row = (r + 1):excess
        if sugarscape(row,c,1) > maxSugar(:)
            maxSugar(m) = sugarscape(row,c,1);
        end
    end
end
% check left
if (c - step) == 0
    if sugarscape(r,column,1) > maxSugar(:)
        maxSugar(m) = sugarscape(r,column,1);
    end
elseif (c - step) < 1
    excess = abs(c - step);
    for column = 1:(c - 1)
        if sugarscape(r,column,1) > maxSugar(:)
            maxSugar(m) = sugarscape(r,column,1);
            m = m + 1;
        end
    end
    for column = (N - excess + 1):N
        if sugarscape(r,column,1) > maxSugar(:)
            maxSugar(m) = sugarscape(r,column,1);
        end
    end
    m = 1;
else
    excess = c - step;
    for column = excess:(c - 1)
        if sugarscape(r,column,1) > maxSugar(:)
            maxSugar(m) = sugarscape(r,column,1);
        end
    end
end
% check right     REMAINS TO BE WRITTEN^^
if (c - step) == 0
    if sugarscape(r,column,1) > maxSugar(:)
        maxSugar(m) = sugarscape(r,column,1);
    end
elseif (c - step) < 1
    excess = abs(c - step);
    for column = 1:(c - 1)
        if sugarscape(r,column,1) > maxSugar(:)
            maxSugar(m) = sugarscape(r,column,1);
            m = m + 1;
        end
    end
    for column = (N - excess + 1):N
        if sugarscape(r,column,1) > maxSugar(:)
            maxSugar(m) = sugarscape(r,column,1);
        end
    end
    m = 1;
else
    excess = c - step;
    for column = excess:(c - 1)
        if sugarscape(r,column,1) > maxSugar(:)
            maxSugar(m) = sugarscape(r,column,1);
        end
    end
end

suga = max(maxSugar);



end
end

s = sugarscape(:,:,1);
imagesc(s);












