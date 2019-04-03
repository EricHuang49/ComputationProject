function [V,ga,gc] = calc_value_policy_w_retirement_mortality(U,agrid,y,ybar,bbeta,r,T,retirement,survivalProb)

na = length(agrid);
ny = length(y.values);

V = zeros(T,na*ny,1);
ga = zeros(T,na*ny,1);
gc = zeros(T,na*ny,1); 

% Age == T
for i = 1 : na*ny
    ia = floor(mod(i-0.05,na))+1;

    asset = agrid(ia);
    income = retirement.ss(T);
    consumption = income + (1+r)*asset;
    V(T,i) = U(consumption);
    if consumption < 0
            V(T,i) = -inf;
    end
    ga(T,i) = min(agrid);
    gc(T,i) = consumption;

%     VV = -inf;
%     for iap = 1 : na
%         ap = agrid(iap);
%         consumption = income + (1+r)*asset - ap;
%         
%         if consumption < 0 || ap < 0
%             break
%         end
%         
%         value = U(consumption) + bbeta*U(ap);
%         if value >= VV
%             VV = value;
%             aChoice = iap;
%             cChoice = consumption;
%         end
%     end
%     V(T,i) = VV;
%     ga(T,i) = aChoice;
%     gc(T,i) = cChoice;
    
end

tic
for age = T-1 : -1 : 1
%     VnextExpected = zeros(na,ny);                   % a and h are T+1 states; r and eepsilon are current 
%     Vnext = V(age+1,:);
%     for iy = 1 : ny
%         Vp = reshape(Vnext,[],ny);
%         Vp = Vp * y.transition(iy,:)';
%         VnextExpected(:,iy) = Vp;
%     end

    Vmatrix = reshape(V(age+1,:),[],ny);
     for i = 1 : na*ny
        ia = floor(mod(i-0.05,na))+1;
        iy = mod(floor((i-0.05)/na),ny)+1;

        asset = agrid(ia);
        if age >= retirement.age
            income = retirement.ss(age);
        else
            income = ybar(age)*y.values(iy);
        end
% 
%         VV = -inf;
%         aChoice = nan;
%         cChoice = nan;
%         for iap = 1 : na
%             ap = agrid(iap);
%             consumption = income + (1+r)*asset - ap;
% 
%             if consumption < 0
%                 break
%             end
% 
%             value = U(consumption) + bbeta*survivalProb(age)*VnextExpected(iap,iy);
% 
%             if value >= VV
%                 VV = value;
%                 aChoice = iap;
%                 cChoice = consumption;
%             end
% 
%         end
% 
%         V(age,i) = VV;
%         ga(age,i) = agrid(aChoice);
%         gc(age,i) = cChoice;

        val = @(a) value_function_w_retirement_mortality(a,Vmatrix,bbeta,U,agrid,y,r,ia,iy,income,survivalProb(age));
        ap = fminbnd(val,min(agrid),max(agrid));
        V(age,i) = -val(ap);
        ga(age,i) = ap;
        gc(age,i) = income + (1+r)*asset - ap;    
    end
    
    finish = toc;
    disp(['Age: ', num2str(age+19), '. Time: ', num2str(finish),' seconds'])

end