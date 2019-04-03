function [aPath,cPath,yPath] = simulate_model(M,Tsim,agrid,y,T,ga,gc)

ny = length(y.values);
na = length(agrid);

aPath = zeros(Tsim+1,M);
cPath = zeros(Tsim,M);

% Generate Markov chains for income
rng(1)
X0 = zeros(1,ny);
%yss_idx = ceil(ny/2);
X0(1) = M;                      % M simulations from state 1
yPath = simulate(dtmc(y.transition),Tsim-1,'X0',X0);

for age = 1 : Tsim    
    currentAsset = aPath(age,:);
    currentIncomeIdx = yPath(age,:);
        
    for j = 1 : M
        if T == inf
            currentAssetPolicy = ga((currentIncomeIdx(j)-1)*na+1:currentIncomeIdx(j)*na);
            currentConsumptionPolicy = gc((currentIncomeIdx(j)-1)*na+1:currentIncomeIdx(j)*na);
        else
            currentAssetPolicy = ga(age,(currentIncomeIdx(j)-1)*na+1:currentIncomeIdx(j)*na);
            currentConsumptionPolicy = gc(age,(currentIncomeIdx(j)-1)*na+1:currentIncomeIdx(j)*na);
        end
    
        aPath(age+1,j) = interpolate(currentAssetPolicy,currentAsset(j),agrid);
        cPath(age,j) = interpolate(currentConsumptionPolicy,currentAsset(j),agrid);
    end
    %cPath(age,:) = y.values(currentIncomeIdx)' + (1+r)*agrid(currentAssetIdx) - agrid(aPath(age+1,:));
end

end