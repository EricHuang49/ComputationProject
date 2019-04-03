clear;clc;
warning off
addpath(genpath('Functions'))
addpath(genpath('../Input'))
outdir = ['..' filesep 'Output' filesep];

%% Set up parameters
bbeta = 1/1.04;
ssigma = 2;
ddelta = 0.95;
ssigmaY = 0.2;
%r = 1/bbeta-1;
r = 0.02;
amin = 0;
amax = 150;

%% Compute the Markov transition matrix for income
ny = 11;
[y.values,y.transition] = mytauchen(0,ddelta,sqrt(1-ddelta^2)*ssigmaY,ny);
y.values = exp(y.values);

%% Specify the grid for assets
na = 300;
agrid = linspace(sqrt(amin),sqrt(amax),na);
agrid = agrid.^2;

%% 2. Value Function Iteration for policy and value functions for infinite periods
% U = @(c) (c.^(1-ssigma)-1)./(1-ssigma);
% T = inf;
% 
% % Calculate value and policy functions
% [V_inf,ga_inf,gc_inf] = calc_value_policy(V0,U,agrid,y,bbeta,r,T);
% 
% % Simulate the model
% M = 1000;                       % number of simulations
% Tsim = 61;
% 
% [aPath_inf,cPath_inf,yPath_inf] = simulate_model(M,Tsim,agrid,y,T,ga_inf,gc_inf);

%% 3. Repeat the exercise for a lifespan of 61 years
% T = 61;
% 
% % Calculate value and policy functions
% [V_finite,ga_finite,gc_finite] = calc_value_policy(V0,U,agrid,y,bbeta,r,T);
% 
% % Simulate the model
% [aPath_finite,cPath_finite,yPath_finite] = simulate_model(M,Tsim,agrid,y,T,ga_finite,gc_finite);

%% 4. Do the exercise for infinite number of periods
T = inf;
ddelta = 0.8;
ssigmaY = [0.2 0.4];
bbeta = 1/1.04;
r = 0.02;

% Specify two income processes with different volatility
[y1.values,y1.transition] = mytauchen(0,ddelta,sqrt(1-ddelta^2)*ssigmaY(1),ny);     % when variance of income innovation is 0.2
y1.values = exp(y1.values);
[y2.values,y2.transition] = mytauchen(0,ddelta,sqrt(1-ddelta^2)*ssigmaY(2),ny);     % when variance of income innovation is 0.4
y2.values = exp(y2.values);

% Log utility
ssigma = 1;
U = @(c) log(c);

% Calculate value and policy functions
V0 = zeros(na*ny,1);
[V4_1,ga4_1,gc4_1] = calc_value_policy(V0,U,agrid,y1,bbeta,r,T);
[V4_2,ga4_2,gc4_2] = calc_value_policy(V0,U,agrid,y2,bbeta,r,T);

% Plot consumption policy functions at average income level
yss_idx = ceil(ny/2);
plot(agrid,gc4_1((yss_idx-1)*na+1:yss_idx*na)');
hold on
plot(agrid,gc4_2((yss_idx-1)*na+1:yss_idx*na)');
legend('g_c, \sigma_y = 0.2','g_c, \sigma_y = 0.4','interpreter','latex')
xlabel('asset')
ylabel('c'' ')
saveas(gcf,[outdir '2_4_consumption.png'])
close(gcf)

% Plot Euler Equation Error
[ee,max_error] = euler_error(bbeta,ssigma,r,y,agrid,ga4_2,gc4_2);
plot(agrid,max_error);
title('Euler Equation Error for VFI');
saveas(gcf,[outdir '2_4_EE.png'])
close(gcf)

% Save the data
save([outdir '2_4.mat']);

%% 5. Repeat question 4 but for finite horizon
T = 61;
r = 0.02;

% Specify two income processes with different volatility
[y1.values,y1.transition] = mytauchen(0,ddelta,sqrt(1-ddelta^2)*ssigmaY(1),ny);     % when variance of income innovation is 0.2
y1.values = exp(y1.values);
[y2.values,y2.transition] = mytauchen(0,ddelta,sqrt(1-ddelta^2)*ssigmaY(2),ny);     % when variance of income innovation is 0.4
y2.values = exp(y2.values);

% Log utility
ssigma = 1;
U = @(c) log(c);

% Calculate value and policy functions
[V5_1,ga5_1,gc5_1] = calc_value_policy(V0,U,agrid,y1,bbeta,r,T);
[V5_2,ga5_2,gc5_2] = calc_value_policy(V0,U,agrid,y2,bbeta,r,T);

% Plot consumption policy functions at ages 0 and 61, at average income level
subplot(1,2,1)
plot(agrid,gc5_1(1,(yss_idx-1)*na+1:yss_idx*na));
hold on
plot(agrid,gc5_2(1,(yss_idx-1)*na+1:yss_idx*na));
title('Age = 20')
legend('g_c, \sigma_y = 0.2','g_c, \sigma_y = 0.4','interpreter','latex')
xlabel('asset')
ylabel('c'' ')
subplot(1,2,2)
plot(agrid,gc5_1(61,(yss_idx-1)*na+1:yss_idx*na));
hold on
plot(agrid,gc5_2(61,(yss_idx-1)*na+1:yss_idx*na));
title('Age = 80')
legend('g_c, \sigma_y = 0.2','g_c, \sigma_y = 0.4','interpreter','latex')
xlabel('asset')
ylabel('c'' ')
saveas(gcf,[outdir '2_5_consumption.png'])
close(gcf)

% Save the data
save([outdir '2_5.mat']);

%% 6.1 Simulate the model and generate paths in life cycles
M = 1000;                       % number of simulations
Tsim = 61;

[aPath1_finite,cPath1_finite,yPath1_finite] = simulate_model(M,Tsim,agrid,y1,T,ga5_1,gc5_1);
cPath1_average = mean(cPath1_finite,2);

[aPath2_finite,cPath2_finite,yPath2_finite] = simulate_model(M,Tsim,agrid,y2,T,ga5_2,gc5_2);
cPath2_average = mean(cPath2_finite,2);

% Plot average consumption paths over the life cycle
plot(20:80,cPath1_average);
hold on
plot(20:80,cPath2_average);
legend('Avg. c, \sigma_y = 0.2','Avg. c, \sigma_y = 0.4','interpreter','latex')
xlabel('age')
ylabel('average consumption')
saveas(gcf,[outdir '2_6_avg_LC_consumption.png'])
close(gcf)

% Save the data
save([outdir '2_6.mat']);

%% 6.2 Experiment with potential ways to introduce hump shape in life cycle consumption

% Set the retirement age to be 70, and set retirement income to be the lowest income level
retirement.age = 50;
retirement.ss = [zeros(retirement.age-1,1); ones(T-retirement.age+1,1)*y.values(1)];
lifeCycleIncome = ones(retirement.age-1,1);
survivalProb = [ones(T-1,1); 0];

% Calculate value and policy functions
[V6,ga6,gc6] = calc_value_policy_w_retirement_mortality(U,agrid,y,lifeCycleIncome,bbeta,r,T,retirement,survivalProb);

% Simulate the model
[aPath6,cPath6,yPath6] = simulate_model(M,Tsim,agrid,y,T,ga6,gc6);

% Plot average life cycle paths for consumption
cPath6_average = mean(cPath6,2);
plot(cPath6_average);
xlabel('age')
ylabel('average consumption')
saveas(gcf,[outdir '2_6_2_avg_LC_consumption_retirement.png'])
close(gcf)

% Save the data
save([outdir '2_6_2.mat']);

%% 7. Introduce retirement and mortality risk
ttheta = 0.7;                       % fraction of average income at age 45 that ages to social security

fileID = fopen('incprofile.txt');
lifeCycleIncome = fscanf(fileID,'%f');
retirement.age = 46;
retirement.ss = [zeros(retirement.age-1,1); ones(T-retirement.age+1,1)*ttheta*lifeCycleIncome(end)];

fileID = fopen('survs.txt');
survivalProb = fscanf(fileID,'%f');

% Calculate value and policy functions
[V7,ga7,gc7] = calc_value_policy_w_retirement_mortality(U,agrid,y,lifeCycleIncome,bbeta,r,T,retirement,survivalProb);

% Simulate the model
[aPath7,cPath7,yPath7] = simulate_model(M,Tsim,agrid,y,T,ga7,gc7);

% Plot average life cycle paths for consumption
cPath7_average = mean(cPath7,2);
plot(cPath7_average);
xlabel('age')
ylabel('average consumption')
saveas(gcf,[outdir '2_7_avg_LC_consumption_retirement.png'])
close(gcf)

% Save the data
save([outdir '2_7.mat'])

%% 8. Improve the fit between the model and the data
% fileID = fopen('consprofile.txt');
% consProfile = textscan(fileID,'%f %f');
% age = consProfile{1};
% cPath_CEX = consProfile{2};
% plot(22:80,cPath_CEX(age<=80 & age==floor(age)));
% hold on
% plot(22:80,cPath7_average(3:end));
% 
% fileID = fopen('hhsize.txt');
% hhSize_5yr = fscanf(fileID,'%f');
% hhAge = [23:5:65 70 78];
% hhSize_lifeCycle = interp1(hhAge,hhSize_5yr,20:1:80,'linear','extrap');
% 
% fileID = fopen('hhequivalence.txt');
% data = textscan(fileID,'%f %f');
% familySize = data{1};
% hhEquiv = data{2};
% 
% hhEquiv_lifeCycle = interp1(familySize,hhEquiv,hhSize_lifeCycle);
% 
% % [aPath8,cPath8,yPath8] = simulate_model(M,Tsim,agrid,y,T,ga8,gc8);
% % cPath8_average = mean(cPath8,2);
% 
% cPath8_hhEquiv = cPath7_average(3:end)./hhEquiv_lifeCycle(3:end)';
% %cPath8_hhEquiv = cPath8_hhEquiv * 1/cPath8_hhEquiv(1);
% plot(22:80,cPath8_hhEquiv);
% hold on
% plot(22:80,cPath_CEX(age<=80 & age==floor(age)));

%% 9. Compute consumption insurance coefficient
ddelta = [0 0.99];
ssigmaY = 0.2;
M = 8000;

% Specify two income processes with different volatility
[y8_1.values,y8_1.transition] = mytauchen(0,ddelta(1),sqrt(1-ddelta(1)^2)*ssigmaY,ny);
y8_1.values = exp(y8_1.values);
[y8_2.values,y8_2.transition] = mytauchen(0,ddelta(2),sqrt(1-ddelta(2)^2)*ssigmaY,ny);
y8_2.values = exp(y8_2.values);

[V8_1,ga8_1,gc8_1] = calc_value_policy(V0,U,agrid,y8_1,bbeta,r,T);
[~,cPath8_1,yPath8_1] = simulate_model(M,Tsim,agrid,y8_1,T,ga8_1,gc8_1);

[V8_2,ga8_2,gc8_2] = calc_value_policy(V0,U,agrid,y8_2,bbeta,r,T);
[~,cPath8_2,yPath8_2] = simulate_model(M,Tsim,agrid,y8_2,T,ga8_2,gc8_2);

% Calculate the consumption insurance coefficient
pphi = nan(T-1,2);
for age = 2 : T
    cov_c_y_1 = cov(log(cPath8_1(age,:))-log(cPath8_1(age-1,:)),log(yPath8_1(age,:))-log(yPath8_1(age-1,:)));
    cov_c_y_2 = cov(log(cPath8_2(age,:))-log(cPath8_2(age-1,:)),log(yPath8_2(age,:))-log(yPath8_2(age-1,:)));
    pphi(age-1,1) = 1 - cov_c_y_1(2)/var(log(yPath8_1(age,:))-log(yPath8_1(age-1,:)));
    pphi(age-1,2) = 1 - cov_c_y_2(2)/var(log(yPath8_2(age,:))-log(yPath8_2(age-1,:)));    
end

plot(pphi(:,1))
hold on
plot(pphi(:,2))
legend('\delta = 0','\delta = 1','interpreter','latex')
xlabel('age')
ylabel('consumption insurance coefficient')
saveas(gcf,[outdir '2_9_consumption_insurance.png'])
close(gcf)

% Save the data
save([outdir '2_9.mat'])

    