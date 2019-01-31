% Piecewise linear interpolation of the invariant distribution
%b = -2; aM = 20; %M  = 200; K=M;
K = 200;
agrid2 = linspace(b,aM,K);
% Endowment process
rho      = 0.2;           % first-order autoregressive coefficient
sigmaint = 0.3;%0.10.3;           % intermediate value to calculate sigma
sigma    = sigmaint*sqrt(1-rho^2); % standard deviation of error_t

% prob is transition matrix of the Markov chain
% logs is the discretized states of log labor earnings
% invdist is the invariant distribution of Markov chain
% alambda and asigma are respectively the theoretical 
% rho and standard deviation of log labor income in Markov chain
N = 2;
[prob,logs,invdist,alambda,asigma]=markovappr(rho,sigma,3,N);
y1 = exp(logs(1));%6
y2 = exp(logs(2));%7
%y = [y1,y2];
% Calculate the inverse of the asset policy function
% the function invaprime @(a)
% Calculate the Invariant distribution of the endowment process
% Choose an initial distribution Lambda0 over the grid of assets a
Lambda0 = zeros(K,N);
for j = 1:N
    for k = 1:K 
        Lambda0(k,j) = (agrid2(k)-b)/(aM-b)*invdist(j);
    end
end
% Make it continous
l0 = @(a,lambda) interp1(agrid2,lambda,a,'linear');
% update distribution for every pair
crit1 = 10^(-10);
dist1 = 1;
lambda = Lambda0;
l1 = zeros(K,N);
a0 = 3;
while dist1>crit1 %dist>crit&&iter<maxiter
        for k = 1:K
            for j = 1:N
                for i = 1:N
                    asol=@(a) fsolve(invaprime(a,Y,r0,a0,cp0,agrid2(k),i),a0);
                    l1(k,j) = l1(k,j)+pi(j,i)*l0(asol,lambda(:,i));
                end
            end
        end

   
    dist1 = max(max(abs(lambda-l1)));
    lambda = l1;
end