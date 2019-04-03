function [V,ga,gc] = calc_value_policy(V0,U,agrid,y,bbeta,r,T,option)

na = length(agrid);
ny = length(y.values);

diff = inf;
tol = 1e-5;

if T == inf
    
    V = V0;
    tempV = zeros(na*ny,1);
    ga = zeros(na*ny,1);
    gc = zeros(na*ny,1);
    
    iter = 1;
    while diff > tol
        
        if exist('option','var') && strcmp(option,'discrete')
            % First calculate the expected continuation value to speed up code
            VnextExpected = zeros(na,ny);                   % a and h are T+1 states; r and eepsilon are current 
            for iy = 1 : ny
                Vp = reshape(V,[],ny);
                Vp = Vp * y.transition(iy,:)';
                VnextExpected(:,iy) = Vp;
            end
        else
            Vmatrix = reshape(V,[],ny);
        end

        if exist('option','var') && strcmp(option,'discrete')
            parfor i = 1 : na*ny
                ia = floor(mod(i-0.05,na))+1;
                iy = mod(floor((i-0.05)/na),ny)+1;

                inc = y.values(iy);
                a = agrid(ia);
                
                VV = -inf;
                aChoice = nan;
                cChoice = nan;

                for iap = 1 : na
                    ap = agrid(iap);
                    c = inc + (1+r)*a - ap;

                    if c < 0
                        break
                    end

                    value = U(c) + bbeta*VnextExpected(iap,iy);

                    if value >= VV
                        VV = value;
                        aChoice = iap;
                        cChoice = c;
                    end

                end

                tempV(i) = VV;
                ga(i) = agrid(aChoice);
                gc(i) = cChoice;
            end
        else
             parfor i = 1 : na*ny
                ia = floor(mod(i-0.05,na))+1;
                iy = mod(floor((i-0.05)/na),ny)+1;

                inc = y.values(iy);
                a = agrid(ia);
                
                val = @(a) value_function(a,Vmatrix,bbeta,U,agrid,y,r,ia,iy);
                ap = fminbnd(val,min(agrid),max(agrid));
                tempV(i) = -val(ap);
                ga(i) = ap;
                gc(i) = inc + (1+r)*a - ap;
            end
        end

        diff = max(abs(V(:)-tempV(:)));
        if (mod(iter,10)==0 || iter ==1)
            fprintf(' Iteration = %d, Sup Diff = %2.8f\n', iter, diff); 
        end
        iter = iter + 1;
        V = tempV;
    end

% Finite-horizon life cycle
else
    
    V = V0;
    ga = zeros(T,na*ny,1);
    gc = zeros(T,na*ny,1); 
    
    for i = 1 : na*ny
        ia = floor(mod(i-0.05,na))+1;
        iy = mod(floor((i-0.05)/na),ny)+1;
        
        asset = agrid(ia);
        income = y.values(iy);
        consumption = income + (1+r)*asset;
        V(T,i) = U(consumption);
        
        if consumption < 0
            V(T,i) = -inf;
        end
        
        ga(T,i) = min(agrid);
        gc(T,i) = consumption;
    end
    
    tic
    for age = T-1 : -1 : 1
%         VnextExpected = zeros(na,ny);                   % a and h are T+1 states; r and eepsilon are current 
%         Vnext = V(age+1,:);
%         for iy = 1 : ny
%             Vp = reshape(Vnext,[],ny);
%             Vp = Vp * y.transition(iy,:)';
%             VnextExpected(:,iy) = Vp;
%         end
        
        Vmatrix = reshape(V(age+1,:),[],ny);
        for i = 1 : na*ny
            ia = floor(mod(i-0.05,na))+1;
            iy = mod(floor((i-0.05)/na),ny)+1;
            
            asset = agrid(ia);
            income = y.values(iy);
            
%             VV = -inf;
%             aChoice = nan;
%             cChoice = nan;
%             for iap = 1 : na
%                 ap = agrid(iap);
%                 consumption = income + (1+r)*asset - ap;
%                 
%                 if consumption < 0
%                     break
%                 end
%                 
%                 value = U(consumption) + bbeta*VnextExpected(iap,iy);
%                 
%                 if value >= VV
%                     VV = value;
%                     aChoice = iap;
%                     cChoice = consumption;
%                 end
%                 
%             end
%            
%             V(age,i) = VV;
%             ga(age,i) = aChoice;
%             gc(age,i) = cChoice;

            val = @(a) value_function(a,Vmatrix,bbeta,U,agrid,y,r,ia,iy);
            ap = fminbnd(val,min(agrid),max(agrid));
            V(age,i) = -val(ap);
            ga(age,i) = ap;
            gc(age,i) = income + (1+r)*asset - ap;
        end
        
        finish = toc;
        disp(['Age: ', num2str(age+19), '. Time: ', num2str(finish),' seconds'])
        
    end
    
end

end