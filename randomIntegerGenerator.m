%% Forest Fires %%
%% Exercise 3.3 %%
%% Kakkos Ioannis %%

%function inverseC = randomIntegerGenerator(alpha,n)
alpha = 1.16;
n = 10^3;
powerLaw = zeros(1,n);
for i = 1:n
    powerLaw(i) = i^(-alpha);
end
C = 1 - powerLaw;
inverseC = C.^(-1);
figure;
histogram(inverseC)
%end


