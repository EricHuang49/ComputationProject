function value = value_function_w_retirement_mortality(ap,Vmatrix,bbeta,U,agrid,y,r,ia,iy,income,survivalProb)
    ny = length(y.values);
    
    % Linearly interpolate the value function
    V_interp = nan(ny,1);
    for i = 1 : ny
        V_interp(i) = interpolate(Vmatrix(:,i),ap,agrid);
    end
    
    c = income + (1+r)*agrid(ia) - ap;
    
    if c < 0
        value = -1e9+c;
    else
        if ~exist('survivalProb','var')
            value = U(c) + bbeta*y.transition(iy,:)*V_interp;
        else
            value = U(c) + bbeta*survivalProb*y.transition(iy,:)*V_interp;
        end
    end
    value = -value;
end