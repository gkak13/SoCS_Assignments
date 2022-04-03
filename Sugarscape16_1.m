%% Simulation of Complex Systems                                %%
%% Kakkos Ioannis                                               %%
%% Exercise 16.1 Schelling's model of segregation               %%
clear all; 

N = 50;
f = 0.1;
nA = 0; nB = 0;
coverage = (1 - f)*N^2;
fulfilled = coverage/2;
grid = zeros(N);
assumptions = false;

%% Create initial grid
while ~ assumptions
    row = randi(N);
    column = randi(N);
    if grid(row,column) == 0
        q = rand;
        if (q < 0.5 && nA < fulfilled)
            grid(row,column) = 1;       % 1 for family A
            nA = nA + 1;
        elseif (q > 0.5 && nB < fulfilled)
            grid(row,column) = 2;       % 2 for family B
            nB = nB + 1;
        end
    end
    if (nA == coverage/2 && nB == coverage/2)
        assumptions = true;
    end
end

%% Store free houses
free = zeros(2,f*N^2);
c = 1;
for i = 1:N
    for j = 1:N
        if grid(i,j) == 0
            free(1,c) = i;
            free(2,c) = j;
            c = c + 1;
        end
    end
end

%% Plotting grid in its initial configuration
figure;
hold on;
title('Initial configuration','Interpreter','Latex');
xlim([1 50]);
ylim([1 50]);
map = [1 1 1; 0 0 1; 0 1 0];
colormap(map);
imagesc(grid);
hold off;

%% Run the rounds
rounds = 1e5;
friends = 0;
happy = zeros(rounds,3);
moveout = zeros(rounds,1);
for i = 1:rounds
    row = randi(N);
    column = randi(N);
    % assume non-periodic boundary conditions
    if grid(row,column) ~= 0    % check that the house is occupied
        % upper left corner
        if (row == N && column == 1)
            if grid(row-1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 1
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % lower left corner
        elseif (row == 1 && column == 1)
            if grid(row+1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 1
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % upper right corner
        elseif (row == N && column == N)
            if grid(row-1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 1
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % lower right corner
        elseif (row == 1 && column == N)
            if grid(row+1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 1
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % lower row   
        elseif row == 1
            if grid(row,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 2
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % left column
        elseif column == 1
            if grid(row+1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 2
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % upper row   
        elseif row == N
            if grid(row,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 2
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % right column   
        elseif column == N
            if grid(row-1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column) == grid(row,column)
                friends = friends + 1;
            end
            if friends < 2
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        % the rest of the grid
        else
            if grid(row-1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column-1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row+1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column+1) == grid(row,column)
                friends = friends + 1;
            end
            if grid(row-1,column) == grid(row,column)
                friends = friends + 1;
            end
            if friends <= 3
                p = randi(length(free));
                x = free(1,p);
                y = free(2,p);
                grid(x,y) = grid(row,column);
                free(1,p) = row;
                free(2,p) = column;
                grid(row,column) = 0;
            end
        end
        friends = 0;
    end
    [happy(i,:),moveout(i)] = happiness(grid);
    moveout(i) = moveout(i)/coverage;
end

%% Fix data for "happy" matrix to account for percentage
happy(:,1) = happy(:,1)./(coverage/2);
happy(:,2) = happy(:,2)./(coverage/2);
happy(:,3) = happy(:,3)./coverage;

%% Plotting grid in its final configuration
figure;
hold on;
title('Final configuration','Interpreter','Latex');
xlim([1 50]);
ylim([1 50]);
map = [1 1 1; 0 0 1; 0 1 0];
colormap(map);
imagesc(grid);
hold off;

%% Plotting happy measure versus relocation events
figure;
hold on;
title('Separate and global happiness versus relocation events','Interpreter','Latex');
xlabel('time','Interpreter','Latex');
ylabel('p','Interpreter','Latex');
ylim([0 1.1]);
plot(happy);
plot(moveout);
hold off;
legend('Happiness A','Happiness B','Global happiness','Relocation events','Interpreter','Latex');

%% Function happiness
function [love,move] = happiness(lattice)
love = zeros(1,3);
N = length(lattice);
friends = 0;
move = 0;
for i = 1:N
    for j = 1:N
        row = i;
        column = j;
        family = lattice(i,j);
        if family ~= 0
        % upper left corner
        if (row == N && column == 1)
            if lattice(row-1,column) == family
                friends = friends + 1;
            end
            if lattice(row,column+1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column+1) == family
                friends = friends + 1;
            end
            if friends >= 1
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % lower left corner    
        elseif (row == 1 && column == 1)
            if lattice(row+1,column) == family
                friends = friends + 1;
            end
            if lattice(row,column+1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column+1) == family
                friends = friends + 1;
            end
            if friends >= 1
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % upper right corner    
        elseif (row == N && column == N)
            if lattice(row-1,column) == family
                friends = friends + 1;
            end
            if lattice(row,column-1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column-1) == family
                friends = friends + 1;
            end
            if friends >= 1
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % lower right corner
        elseif (row == 1 && column == N)
            if lattice(row+1,column) == family
                friends = friends + 1;
            end
            if lattice(row,column-1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column-1) == family
                friends = friends + 1;
            end
            if friends >= 1
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % lower row   
        elseif row == 1
            if lattice(row,column-1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column-1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column) == family
                friends = friends + 1;
            end
            if lattice(row+1,column+1) == family
                friends = friends + 1;
            end
            if lattice(row,column+1) == family
                friends = friends + 1;
            end
            if friends >= 2
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % left column
        elseif column == 1
            if lattice(row+1,column) == family
                friends = friends + 1;
            end
            if lattice(row+1,column+1) == family
                friends = friends + 1;
            end
            if lattice(row,column+1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column+1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column) == family
                friends = friends + 1;
            end
            if friends >= 2
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % upper row   
        elseif row == N
            if lattice(row,column-1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column-1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column) == family
                friends = friends + 1;
            end
            if lattice(row-1,column+1) == family
                friends = friends + 1;
            end
            if lattice(row,column+1) == family
                friends = friends + 1;
            end
            if friends >= 2
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % right column   
        elseif column == N
            if lattice(row-1,column) == family
                friends = friends + 1;
            end
            if lattice(row-1,column-1) == family
                friends = friends + 1;
            end
            if lattice(row,column-1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column-1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column) == family
                friends = friends + 1;
            end
            if friends >= 2
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        % the rest of the grid
        else
            if lattice(row-1,column-1) == family
                friends = friends + 1;
            end
            if lattice(row,column-1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column-1) == family
                friends = friends + 1;
            end
            if lattice(row+1,column) == family
                friends = friends + 1;
            end
            if lattice(row+1,column+1) == family
                friends = friends + 1;
            end
            if lattice(row,column+1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column+1) == family
                friends = friends + 1;
            end
            if lattice(row-1,column) == family
                friends = friends + 1;
            end
            if friends > 3
                love(1,family) = love(1,family) + 1;
            else
                move = move + 1;
            end
        end
        end
        friends = 0;
    end
end
love(1,3) = love(1,1) + love(1,2);
end
