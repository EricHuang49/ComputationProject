clear;clc;
warning off
addpath(genpath('Functions'))
addpath(genpath('../Input'))
outdir = ['..' filesep 'Output' filesep];

%% Set up parameters
rrho = 0.04;
bbeta = 1/(1+rrho);
ssigma = 1.5;
ddelta = 0.08;              % depreciation rate
aalpha = 0.36;              % capital share in Cobb-Douglas
ttheta = 0.9;               % coefficient for labor endowment process
ssigma_epsilon = 0.4;       % stdev of labor endowment innovations
amin = 0;
amax = 50;

U = @(c) (c.^(1-ssigma)-1)./(1-ssigma);

%% Compute Markov Transition matrix for labor efficiency
ny = 7;
[y.values,y.transition] = mytauchen(0,ttheta,sqrt(1-ttheta^2)*ssigma_epsilon,ny);
y.values = exp(y.values);

% Find out invariant distribution for labor endowment
[V,~] = eig(y.transition');
pi_y = abs(V(:,1)) / sum(abs(V(:,1)));

% Average labor endowment
y_avg = pi_y' * y.values;

% Normalize labor endowment
y.values = y.values / y_avg;

%% Generate grid for asset
na = 500;
% agrid = linspace(sqrt(amin),sqrt(amax),na);
% agrid = agrid.^2;
agrid = linspace(amin,amax,na);
ttau = 0;
llambda = 0;

f1 = @(r) excess_capital(r,aalpha,bbeta,ddelta,U,y,agrid,ttau,llambda);
r1 = fminbnd(f1,-ddelta,rrho);

%% Add tax system
ttau = 0.18;
llambda = 0.11;

f2 = @(r) excess_capital(r,aalpha,bbeta,ddelta,U,y,agrid,ttau,llambda);
r2 = fminbnd(f2,-ddelta,rrho);

tax_revenue = calc_tax(r,aalpha,ddelta,y,agrid,ttau,llambda);