function value = value_function(ap,Vmatrix,bbeta,U,agrid,y,r,ia,iy)
    ny = length(y.values);
    
    % Linearly interpolate the value function
    V_interp = nan(ny,1);
    for i = 1 : ny
        V_interp(i) = interpolate(Vmatrix(:,i),ap,agrid);
    end
    
    c = y.values(iy) + (1+r)*agrid(ia) - ap;
    
    if c < 0
        value = -1e9+c;
    else   
        value = U(c) + bbeta*y.transition(iy,:)*V_interp;
    end
    value = -value;
end