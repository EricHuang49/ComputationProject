function tax = calc_tax(r,aalpha,ddelta,y,agrid,ttau,llambda)
    na = length(agrid);
    ny = length(y.values);

    [W,~] = eig(y.transition');
    pphi = abs(W(:,1))/sum(abs(W(:,1)));
    
    capital_demand = ((r+ddelta)/aalpha).^(1/(aalpha-1));
    wage = (1-aalpha) * capital_demand^aalpha;
    
    yvec = wage * ones(na,1) * y.values';
    yvec = reshape(yvec,na*ny,1);
    after_tax_yvec = (1-ttau).*yvec.^(1-llambda);
    tax = pphi' * (yvec-after_tax_yvec);
end