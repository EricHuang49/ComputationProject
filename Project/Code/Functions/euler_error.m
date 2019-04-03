function [euler_error,max_error] = euler_error(bbeta,ssigma,r,y,agrid,ga,gc)

ny = length(y.values);
na = length(agrid);

iy_mid = ceil(ny/2);
euler_error = zeros(1,na);
       
for ia = 1 : na
    c = gc((iy_mid-1)*na+ia);
    ap = ga((iy_mid-1)*na+ia);
    
    cp = zeros(ny,1);
    for iyp = 1 : ny
        c_values = gc((iyp-1)*na+1:iyp*na);
        cp(iyp) = interpolate(c_values,ap,agrid);
    end
    
    euler_lhs = bbeta*(1+r)*(y.transition(iy_mid,:)*(cp./c).^(-ssigma));
    euler_rhs = 1;
    euler_error(ia) = 1 - euler_rhs/euler_lhs;
end
euler_error = log10(abs(euler_error));
max_error = max(euler_error,[],1);