function Q = transition_matrix(ga,y,agrid)
    na = length(agrid);
    ny = length(y.values);
    
    Q = zeros(na*ny,na*ny);
    
    parfor i = 1 : na*ny       
        iy = mod(floor((i-0.05)/na),ny)+1;
        
        ap = ga(i);
        adiff = abs(agrid-ap);
        
        temp = (adiff==min(adiff))'*y.transition(iy,:);
        Q(i,:) = reshape(temp,1,na*ny);
    end
    
end