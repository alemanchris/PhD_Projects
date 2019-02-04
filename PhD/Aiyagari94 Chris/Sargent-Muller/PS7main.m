%%% this solves ECON 5300, 2015, Ex 7.1(f)-(h)

%% 0. housekeeping
clear all;
clc;

%% 1. parameters and functional forms

% parameters
alpha  = 0.36;%1/3;                   % capital income share
b      = 0;                     % borrowing constraint
beta   = 0.96;%97/100;                % subjective discount factor
delta  = 0.08;%5/100;                 % depreciation rate of physical capital
gamma  = 3;%3/2;                   % inverse elasticity of intertemporal substitution
varphi = 0.2;%%2/3                   % Frisch elasticity of labor supply
%rho    = 5/10;                  % persistence parameter prodictivity
N      = 2;                     % number of possible productivity realizations
%y1     = 95/100;                % low productivity realization
%y2     = 105/100;               % high productivity realization

% transition probability matrix (productivity)
%pi      = zeros(N,N);

% probabilities have to sum up to 1
%{
pi(1,1) = rho;                  % transition from state 1 to state 1
pi(1,2) = 1-rho;                % transition from state 2 to state 1
pi(2,1) = 1-rho;                % transition from state 1 to state 2
pi(2,2) = rho;                  % transition from state 2 to state 2
%}
% (inverse) marginal utility functions
up    = @(c) c.^(-gamma);        % marginal utility of consumption
invup = @(x) x.^(-1/gamma);      % inverse of marginal utility of consumption
%vp    = @(h) h.^(1/varphi);      % marginal disutility of labor supply
invvp = @(x) x.^(varphi);        % inverse marginal disutility of labor supply   
%% 1.2 Income process Sargent

rho_auto = 0.2;           % first-order autoregressive coefficient
sigmaint = 0.4;           % intermediate value to calculate sigma
sigma    = sigmaint*sqrt(1-rho_auto^2); % standard deviation of error_t

% prob is transition matrix of the Markov chain
% logs is the discretized states of log labor earnings
% invdist is the invariant distribution of Markov chain
% alambda and asigma are respectively the theoretical 
% rho and standard deviation of log labor income in Markov chain
%
[pi,logs,invdist2,alambda,asigma]=markovappr(rho_auto,sigma,3,N);
%pi=pi';
s_income = exp(logs);
y1 = s_income(1);
y2 = s_income(2);
%}
%% 2. discretization

% set up asset grid
M  = 250;%250;                        % number of asset grid points
aM = 50;  %50                        % maximum asset level
A  = linspace(b,aM,M)';          % equally-spaced asset grid from a_1=b to a_M

% set up productivity grid
Y  = [y1,y2]';                   % grid for productivity

% vectorize the grid in two dimensions
Amat = repmat(A,1,N);            % values of A change vertically
Ymat = repmat(Y',M,1);           % values of Y change horizontally

% % this is the built-in alternative
% [Amat,Ymat] = ndgrid(A,Y);

%% 3. endogenous functions

% optimal labor supply
H  = @(c,y,w) invvp(up(c).*y.*w);

% current consumption level, cp0(anext,ynext) is the guess
C0 = @(cp0,r) invup(beta*(1+r)*up(cp0)*pi');
                
% current asset level, c0 = C0(cp0(anext,ynext))
A0 = @(anext,y,c0,r,w) 1/(1+r)...
                    *(c0+anext-H(c0,y,w).*y.*w);
                
%% 4. solve for the stationary equilibrium

% convergence criterion for consumption iteration
crit = 10^(-6);

% parameters of the simulation
% note: 
% it takes quite some simulation periods to get to the stationary
% distribution. Choose a high T >= 10^(-4) once the algorithm is running.
I = 10^(4);             % number of individuals
T = 10^(4);             % number of periods

% choose interval where to search for the stationary interest rate
% note: 
% the staionary distribution is very sensitive to the interst rate. 
% make use of the theoretical result that the stationary rate is slightly 
% below 1/beta-1
r0  = (1/beta-1)-[10^(-12),10^(-4)];
r0 = [0.001, 0.316];
% set up an anonymous function
fprintf('Start solving the Aiyagari model... \n');

tic;
myfun   = @(r) stationary_equilibrium(r,crit,I,T,Amat,Ymat,alpha,b,delta,pi(1,1),varphi,A0,C0,H);
%{
options = optimset('display','iter','TolX',1e-8,'MaxIter',20);
rstar   = fzero(myfun,r0,options);
fprintf('Done with the Aiyagari model in %f sec. \n',toc);
%}
% get the simulated asset levels
fprintf('Fetching the wealth distribution... \n');
rstar = 0.036;
[r,at] = stationary_equilibrium(rstar,crit,I,T,Amat,Ymat,alpha,b,delta,pi(1,1),varphi,A0,C0,H);

%% 5. plot the wealth distribution

% use the last 100 periods
[n,xout] = hist(at(:,T-100:T),M);     % choose M bins
bar(xout,n/sum(n));                   % relative frequency is n/sum(n)
title(sprintf('Wealth distribution, rho = %2.1f',pi(1,1)));
xlabel('Asset level')
ylabel('Relative frequency')
