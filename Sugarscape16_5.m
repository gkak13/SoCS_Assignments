%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis                                               %%
%% Exercise 16.5 A simple Sugarscape                            %%
clear all;

N = 50;
A = 400;
agents = zeros(3,A);
rounds = 5e2;

%% Initialize vision, metabolic rate and initial sugar endowment for each 
%% agent
for i = 1:A
    agents(1,i) = unidrnd(6);
    agents(2,i) = unidrnd(4);
    agents(3,i) = unidrnd(21) + 4;
end
visionInitial = agents(1,:);
vision = agents(1,:);
metabolismInitial = agents(2,:);
metabolism = agents(2,:);
sugar = zeros(rounds,A);
sugar(1,:) = agents(3,:);

%% Initialize sugarscape
sugarscape = zeros(N,N,2);
x01 = 13; y01 = 36; x02 = 36; y02 = 13;
r1 = 5; r2 = 4; r3 = 4; r4 = 4;
t1 = 1; t2 = 1; t3 = 1; t4 = 1;
for i = 1:N
for j = 1:N
    dif1 = sqrt((i - x01)^2 + (j - y01)^2);
    dif2 = sqrt((i - x02)^2 + (j - y02)^2);
    if dif1 <= r1
        sugarscape(i,j) = 4;
        regenerationCells1(t1,:) = [i j];
        t1 = t1 + 1;
    elseif (dif1 > r1 && dif1 <= (r1 + r2))
        sugarscape(i,j) = 3;
        regenerationCells2(t2,:) = [i j];
        t2 = t2 + 1;
    elseif (dif1 > (r1 + r2) && dif1 <= (r1 + r2 + r3))
        sugarscape(i,j) = 2;
        regenerationCells3(t3,:) = [i j];
        t3 = t3 + 1;
    elseif (dif1 > (r1 + r2 + r3) && dif1 <= (r1 + r2 + r3 + r4))
        sugarscape(i,j) = 1;
        regenerationCells4(t4,:) = [i j];
        t4 = t4 + 1;
    end
    if dif2 <= r1
        sugarscape(i,j) = 4;
        regenerationCells1(t1,:) = [i j];
        t1 = t1 + 1;
    elseif (dif2 > r1 && dif2 <= (r1 + r2))
        sugarscape(i,j) = 3;
        regenerationCells2(t2,:) = [i j];
        t2 = t2 + 1;
    elseif (dif2 > (r1 + r2) && dif2 <= (r1 + r2 + r3))
        sugarscape(i,j) = 2;
        regenerationCells3(t3,:) = [i j];
        t3 = t3 + 1;
    elseif (dif2 > (r1 + r2 + r3) && dif2 <= (r1 + r2 + r3 + r4))
        sugarscape(i,j) = 1;
        regenerationCells4(t4,:) = [i j];
        t4 = t4 + 1;
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
chosen = zeros(1,A);
for i = 1:(rounds - 1)
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
maxSugar = 0;
% check upwards
if (r - step) == 0
    if (sugarscape(N,c,1) > maxSugar && sugarscape(N,c,2) == 0)
        maxSugar = sugarscape(N,c,1);
        maxSugar_r = N;
        maxSugar_c = c;
    end
    if r ~= 1
        for row = 1:(r - 1)
            if (sugarscape(row,c,1) > maxSugar && sugarscape(row,c,2) == 0)
                maxSugar = sugarscape(row,c,1);
                maxSugar_r = row;
                maxSugar_c = c;
            end
        end
    end
elseif (r - step) < 0
    excess = abs(r - step);
    for row = (N - excess):N
        if (sugarscape(row,c,1) > maxSugar && sugarscape(row,c,2) == 0)
            maxSugar = sugarscape(row,c,1);
            maxSugar_r = row;
            maxSugar_c = c;
        end
    end
    for row = 1:(r - 1)
        if (sugarscape(row,c,1) > maxSugar && sugarscape(row,c,2) == 0)
            maxSugar = sugarscape(row,c,1);
            maxSugar_r = row;
            maxSugar_c = c;
        end
    end
else
    excess = r - step;
    for row = excess:(r - 1)
        if (sugarscape(row,c,1) > maxSugar && sugarscape(row,c,2) == 0)
            maxSugar = sugarscape(row,c,1);
            maxSugar_r = row;
            maxSugar_c = c;
        end
    end
end
% check downwards
if (r + step) > N
    excess = r + step - N;
    for row = 1:excess
        if (sugarscape(row,c,1) > maxSugar && sugarscape(row,c,2) == 0)
            maxSugar = sugarscape(row,c,1);
            maxSugar_r = row;
            maxSugar_c = c;
        end
    end
    for row = (r + 1):N
        if (sugarscape(row,c,1) > maxSugar && sugarscape(row,c,2) == 0)
            maxSugar = sugarscape(row,c,1);
            maxSugar_r = row;
            maxSugar_c = c;
        end
    end
else
    excess = r + step;
    for row = (r + 1):excess
        if (sugarscape(row,c,1) > maxSugar && sugarscape(row,c,2) == 0)
            maxSugar = sugarscape(row,c,1);
            maxSugar_r = row;
            maxSugar_c = c;
        end
    end
end
% check left
if (c - step) == 0
    if (sugarscape(r,N,1) > maxSugar && sugarscape(r,N,2) == 0)
        maxSugar = sugarscape(r,N,1);
        maxSugar_r = r;
        maxSugar_c = N;
    end
    for column = 1:(c - 1)
        if (sugarscape(r,column,1) > maxSugar && sugarscape(r,column,2) == 0)
            maxSugar = sugarscape(r,column,1);
            maxSugar_r = r;
            maxSugar_c = column;
        end
    end
elseif (c - step) < 0
    excess = abs(c - step);
    for column = 1:(c - 1)
        if (sugarscape(r,column,1) > maxSugar && sugarscape(r,column,2) == 0)
            maxSugar = sugarscape(r,column,1);
            maxSugar_r = r;
            maxSugar_c = column;
        end
    end
    for column = (N - excess):N
        if (sugarscape(r,column,1) > maxSugar && sugarscape(r,column,2) == 0)
            maxSugar = sugarscape(r,column,1);
            maxSugar_r = r;
            maxSugar_c = column;
        end
    end
else
    excess = c - step;
    for column = excess:(c - 1)
        if (sugarscape(r,column,1) > maxSugar && sugarscape(r,column,2) == 0)
            maxSugar = sugarscape(r,column,1);
            maxSugar_r = r;
            maxSugar_c = column;
        end
    end
end
% check right  
if (c + step) <= N
    for column = (c + 1):(c + step)
        if (sugarscape(r,column,1) > maxSugar && sugarscape(r,column,2) == 0)
            maxSugar = sugarscape(r,column,1);
            maxSugar_r = r;
            maxSugar_c = column;
        end
    end
else
    excess = c + step - N;
    for column = 1:excess
        if (sugarscape(r,column,1) > maxSugar && sugarscape(r,column,2) == 0)
            maxSugar = sugarscape(r,column,1);
            maxSugar_r = r;
            maxSugar_c = column;
        end
    end
    for column = (c + 1):N
        if (sugarscape(r,column,1) > maxSugar && sugarscape(r,column,2) == 0)
            maxSugar = sugarscape(r,column,1);
            maxSugar_r = r;
            maxSugar_c = column;
        end
    end
end

if maxSugar > 0  % check that high sugar content was indeed found!
    
    % Agent moves to higher sugar content
    sugarscape(maxSugar_r,maxSugar_c,2) = agent;
    % Initial position of agent emptied
    sugarscape(r,c,2) = 0;
    
    % Perform sugar content variation function
    sugar(i + 1,agent) = sugar(i,agent) +...
        sugarscape(maxSugar_r,maxSugar_c,1) - metabolism(agent);
    
    % Check the agent's sugar load
    if sugar(i + 1,agent) <= 0
        sugarscape(maxSugar_r,maxSugar_c,2) = 0;  % kill the agent at
                                                  % the new cell it arrived
        metabolism(agent) = 0;
        vision(agent) = 0;
        sugar(i + 1,agent) = 0;
    end
    
    % Sugar is harvested from the cell the agent visited
    sugarscape(maxSugar_r,maxSugar_c,1) = 0;
    
end

end
%% Sugar regrowth rule
%  First layer, 4 units of sugar max
for k = 1:length(regenerationCells1)
    x = regenerationCells1(k,1);
    y = regenerationCells1(k,2);
    if sugarscape(x,y,1) <= 3
        sugarscape(x,y,1) = sugarscape(x,y,1) + 1;
    end
end
%  Second layer, 3 units of sugar max
for k = 1:length(regenerationCells2)
    x = regenerationCells2(k,1);
    y = regenerationCells2(k,2);
    if sugarscape(x,y,1) <= 2
        sugarscape(x,y,1) = sugarscape(x,y,1) + 1;
    end
end
%  Third layer, 2 units of sugar max
for k = 1:length(regenerationCells3)
    x = regenerationCells3(k,1);
    y = regenerationCells3(k,2);
    if sugarscape(x,y,1) <= 1
        sugarscape(x,y,1) = sugarscape(x,y,1) + 1;
    end
end
%  Fourth layer, 1 unit of sugar max
for k = 1:length(regenerationCells4)
    x = regenerationCells4(k,1);
    y = regenerationCells4(k,2);
    if sugarscape(x,y,1) == 0   
        sugarscape(x,y,1) = sugarscape(x,y,1) + 1;
    end
end

end

% Calculate the number of alive agents
m = 0;
for i = 1:N
    for j = 1:N
        if sugarscape(i,j,2) ~= 0
        m = m + 1;
        end
    end
end
X = ['The number of alive agents is: ',num2str(m)];
disp(X)

s1 = sugarscape(:,:,1);
s2 = sugarscape(:,:,2);

sug1 = sugar(1,:);
sug21 = sugar(21,:);
sug41 = sugar(41,:);
sug61 = sugar(61,:);

figure;
hold on;
title('Sugar content','Interpreter','Latex');
xlabel('x - coordinates of Lattice','Interpreter','Latex');
ylabel('y - coordinates of Lattice','Interpreter','Latex');
xlim([1 50]);
ylim([1 50]);
map = [1 1 1; .929 .694 .125; .85 .325 .098; .635 .078 .184; 1 0 0];
colormap(map);
imagesc(s1);
hold off;

figure;
hold on;
title('Agents in sugarscape','Interpreter','Latex');
xlabel('x - coordinates of Lattice','Interpreter','Latex');
ylabel('y - coordinates of Lattice','Interpreter','Latex');
xlim([1 50]);
ylim([1 50]);
map = [1 1 1; .466 .674 .188];
colormap(map);
imagesc(s2);
hold off;

figure;
hold on;
title('Metabolism','Interpreter','Latex');
xlabel('m (metabolism)','Interpreter','Latex');
ylabel('N','Interpreter','Latex');
histogram(metabolismInitial);
histogram(metabolism,'BinLimits',[.5 4.5]);
hold off;
legend('Initial metabolism','Final metabolism','Interpreter','Latex');

figure;
hold on;
title('Vision','Interpreter','Latex');
xlabel('v (vision/velocity)','Interpreter','Latex');
ylabel('N','Interpreter','Latex');
histogram(visionInitial);
histogram(vision,'BinLimits',[.5 6.5]);
hold off;
legend('Initial vision','Final vision','Interpreter','Latex');

figure;
hold on; 
title('Wealth distribution','Interpreter','Latex');
xlabel('s (sugar/wealth)','Interpreter','Latex');
ylabel('N');
histogram(sug1,5,'FaceColor',[0 0.4470 0.7410]);
histogram(sug41,'BinLimits',[.5 160],'FaceColor',[0.9290 0.6940 0.1250]);
hold off;
legend('Round 1','Round 41','Interpreter','Latex')

figure;
hold on; 
title('Wealth distribution','Interpreter','Latex');
xlabel('s (sugar/wealth)','Interpreter','Latex');
ylabel('N');
histogram(sug21,'BinLimits',[.5 160],'FaceColor',[0.4660 0.6740 0.1880]);
histogram(sug61,'BinLimits',[.5 160],'FaceColor',[0.8500 0.3250 0.0980]);
hold off;
legend('Round 21','Round 61','Interpreter','Latex')

figure;
hold on; 
title('Wealth distribution','Interpreter','Latex');
xlabel('s (sugar/wealth)','Interpreter','Latex');
ylabel('N');
histogram(sug1,5,'FaceColor',[0 0.4470 0.7410]);
histogram(sug21,'BinLimits',[.5 160],'FaceColor',[0.4660 0.6740 0.1880]);
histogram(sug41,'BinLimits',[.5 160],'FaceColor',[0.9290 0.6940 0.1250]);
histogram(sug61,'BinLimits',[.5 160],'FaceColor',[0.8500 0.3250 0.0980]);
hold off;
legend('Round 1','Round 21','Round 41','Round 61','Interpreter','Latex')
