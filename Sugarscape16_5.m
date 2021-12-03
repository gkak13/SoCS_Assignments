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

sugarscape = zeros(N,N,2);
for i = 1:N/2
    























