function [V2,ga,gc] = calc_value_policy_GE(r,aalpha,bbeta,ddelta,U,y,agrid,ttau,llambda)

    na = length(agrid);
    ny = length(y.values);
    
    capital_demand = ((r+ddelta)/aalpha).^(1/(aalpha-1));
    wage = (1-aalpha) * capital_demand^aalpha;
        
    after_tax_income.values = (1-ttau).*(y.values .* wage).^(1-llambda);
    after_tax_income.transition = y.transition;
    
    T = inf;
    
    % Use VFI with grid search (100 asset grid points) to get an approximate value function
    agrid_small = linspace(min(agrid),max(agrid),100);
    V0 = zeros(100*ny,1);  
    [V,~,~] = calc_value_policy(V0,U,agrid_small,after_tax_income,bbeta,r,T,'discrete');
    
    % Interpolate the value function on the finner grid
    Vmatrix = reshape(V,[],ny);
    V1 = nan(na*ny,1);
    for ia = 1 : na
        for iy = 1 : ny
            idx = (iy-1)*na+ia;
            asset = agrid(ia);
            V1(idx) = interpolate(Vmatrix(:,iy),asset,agrid_small);
        end
    end
    
    % Use VFI with linear interpolation to get a precise policy function
    [V2,ga,gc] = calc_value_policy(V1,U,agrid,y,bbeta,r,T);
end