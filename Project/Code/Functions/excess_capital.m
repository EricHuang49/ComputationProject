function d = excess_capital(r,aalpha,bbeta,ddelta,U,y,agrid,ttau,llambda)
    fprintf(' r = %2.8f\n', r); 
    
    na = length(agrid);
    ny = length(y.values);
    
    [~,ga,~] = calc_value_policy_GE(r,aalpha,bbeta,ddelta,U,y,agrid,ttau,llambda);
    
     % Cast the solution into a Markov chain
    Q = transition_matrix(ga,y,agrid);
    
    % Find the invariant distribution    
%     pphi = ones(na*ny,1)/(na*ny);
%     diff = inf;
%     tol = 1e-6;
%     
%     while diff > tol
%         pphi_new = Q' * pphi;
%         diff = sum(abs(pphi-pphi_new));
%         pphi = pphi_new;
%     end

    [W,~] = eig(Q');
    pphi = abs(W(:,1))/sum(abs(W(:,1)));
    
    avec = agrid' * ones(1,ny);
    avec = reshape(avec,na*ny,1);
    capital_supply = pphi' * avec;
    
    capital_demand = ((r+ddelta)/aalpha).^(1/(aalpha-1));
    
    d = abs(capital_supply - capital_demand);
    fprintf(' Excess asset Supply = %2.8f\n', capital_supply - capital_demand);
end