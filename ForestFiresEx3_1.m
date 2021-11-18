%% Forest Fires %%
%% Exercise 3.1 %%
%% Kakkos Ioannis %%

%clear all; clc;

N = 16;
p = 0.01;                 % tree growth parameter
f = 0.2;                  % lightning strike probability
timesteps = 10^3;
grid = zeros(N);
fires = 0;
treesBorn = 0;
k = 0;
fireStarter = zeros(1,timesteps);

for e = 1:timesteps

    % fill the grid with trees with probability p (typewriter scheme)
    for i = 1:N
        for j = 1:N
            r = rand;
            if (r <= p && grid(i,j) == 0)
                treesBorn = treesBorn + 1;
                grid(i,j) = 1;
            end
        end
    end

    % simulate the lightning strike probability
    q = rand;
    row = randi(N);
    column = randi(N);
    if (q <= f && grid(row,column) == 1)
        fires = fires + 1;
        grid(row,column) = 2;
        k = k + 1;
        % simulate the fire spread
        % record the fire sizes for each time a tree was struck
        fireStarter(k) = fireStarter(k) + 1;
        for i = 2:N-1
            for j = 2:N-1
                if grid(i,j) == 2
                    if grid(i+1,j) == 1
                        grid(i+1,j) = 2;
                        fireStarter(k) = fireStarter(k) + 1;
                    end
                    if grid(i-1,j) == 1
                        grid(i-1,j) = 2;
                        fireStarter(k) = fireStarter(k) + 1;
                    end
                    if grid(i,j+1) == 1
                        grid(i,j+1) = 2;
                        fireStarter(k) = fireStarter(k) + 1;
                    end
                    if grid(i,j-1) == 1
                        grid(i,j-1) = 2;
                        fireStarter(k) = fireStarter(k) + 1;
                    end
                end
            end
        end
        % simulate the boundary conditions
        % upper and lower row
        for j = 1:N
            if grid(1,j) == 2
                grid(N,j) = 2;
                fireStarter(k) = fireStarter(k) + 1;
            elseif grid(N,j) == 2
                grid(1,j) = 2;
                fireStarter(k) = fireStarter(k) + 1;
            end
        end
        % right and left column
        for i = 1:N
            if grid(i,1) == 2
                grid(i,N) = 2;
                fireStarter(k) = fireStarter(k) + 1;
            elseif grid(i,N) == 2
                grid(i,1) = 2;
                fireStarter(k) = fireStarter(k) + 1;
            end
        end
    end

    % set the burnt cells equal to zero
    grid = resetGrid(grid,N);

    % plot the results after zeroing the burnt cells
%     if rem(e,100) == 1
%         figure;
%         pcolor(grid);
%         colormap(summer(2));
%     end
end

% figure;
% histogram(fireStarter);
% hold on;
% title('Fire size distribution histogram');
% xlabel('T');
% ylabel('Fire size');
% hold off;

sumOfTrees = sum(fireStarter);


function x = resetGrid(array,n)
x = array;
for i = 1:n
    for j = 1:n
        if x(i,j) == 2
            x(i,j) = 0;
        end
    end
end
end



