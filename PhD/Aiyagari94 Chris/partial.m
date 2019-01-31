%function [arule]=partial(r0,Amat,Ymat,alpha,b,delta,rho,varphi)
%stationary_equilibrium(r0,crit,I,T,Amat,Ymat,alpha,b,delta,rho,varphi)
clear all
clc
close all
% convergence criterion for consumption iteration
crit = 10^(-6);
% maximum iterations (avoids infinite loops)
maxiter = 10^(3);%10^(3);
%% 1. parameters and functional forms
r0 = 0.015; %0.04 0.015
% parameters
alpha  = 0.36;                   % capital income share
b      = -36;  %OR-2 -35                  % borrowing constraint
beta   = 0.97;%97/100;                % subjective discount factor
delta  = 0.08;%5/100;                 % depreciation rate of physical capital
gamma  = 3;%mu%1.13/2;              % inverse elasticity of intertemporal substitution
%varphi = 0.8;%1/3;%0.5;%2/3;                   % Frisch elasticity of labor supply
%rho    = 5/10;                  % persistence parameter prodictivity
N      = 2;                     % number of possible productivity realizations

rho      = 0.2;           % first-order autoregressive coefficient
sigmaint = 0.3;%0.10.3;           % intermediate value to calculate sigma
sigma    = sigmaint*sqrt(1-rho^2); % standard deviation of error_t

% prob is transition matrix of the Markov chain
% logs is the discretized states of log labor earnings
% invdist is the invariant distribution of Markov chain
% alambda and asigma are respectively the theoretical 
% rho and standard deviation of log labor income in Markov chain

[prob,logs,invdist,alambda,asigma]=markovappr(rho,sigma,3,N);
y1 = exp(logs(1));%6
y2 = exp(logs(2));%7


%y1     = 95/100;%95/100;                % low productivity realization
%y2     = 105/100;%105/100;               % high productivity realization

%transition probability matrix (productivity)
pi      = zeros(N,N);

% probabilities have to sum up to 1
pi(1,1) = prob(1,1);%rho;                  % transition from state 1 to state 1
pi(1,2) = prob(1,2);%1-rho;                % transition from state 2 to state 1
pi(2,1) = prob(2,1);%1-rho;                % transition from state 1 to state 2
pi(2,2) = prob(2,2);%rho;                  % transition from state 2 to state 2

% (inverse) marginal utility functions
up    = @(c) c.^(-gamma);        % marginal utility of consumption
invup = @(x) x.^(-1/gamma);      % inverse of marginal utility of consumption
%vp    = @(h) h.^(1/varphi);      % marginal disutility of labor supply
%invvp = @(x) x.^(varphi);        % inverse marginal disutility of labor supply    

%% 2. discretization

% set up asset grid
M  = 150;%150250 or 400 good;                        % number of asset grid points
aM = 70;%20;%45;                         % maximum asset level
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
H  = 1;%@(c,y,w) invvp(up(c).*y.*w);

% current consumption level, cp0(anext,ynext) is the guess
C0 = @(cp0,r) invup(beta*(1+r)*up(cp0)*pi');
                
% current asset level, c0 = C0(cp0(anext,ynext))
A0 = @(anext,y,c0,r,w) 1/(1+r)...
                    *(c0+anext-H.*y.*w);
                
%% 4. solve for the stationary equilibrium


%%
% this function
% (1) solves for the consumption decision rule, given an
% intereste rate, r0
% (2) simulates the stationary equilibrium associated with the interest
% rate, r0
% (3) (i) returns the residual between the given interes rate r0 and the one
% implied by the stationary aggregate capital and labor supply.
% (ii) returns as an optional output the wealth distribution.

% get dimensions of the grid
[M,N] = size(Amat);

% get productivity realizations from the first row
y1 = Ymat(1,1);
y2 = Ymat(1,2);

% compute the wage according to marginal pricing and a Cobb-Douglas 
% production function
w0 = (1-alpha)*(alpha/(r0+delta))^(alpha/(1-alpha));
% stationary rate is slightly below 1/beta-1?
if r0>=1/beta-1
    disp('to high interest rate')
elseif w0*y1/r0<b
    disp('Violates natural borrowing limit')
else
    disp('Alles gut')
end
%w = (1-alpha)*(alpha/(r0+delta))^alpha)^(1/(1-alpha);
% initial guess on future consumption (consume asset income plus labor
% income from working h=1.
cp0 = r0*Amat+Ymat*w0;

%%% iterate on the consumption decision rule

% distance between two successive functions
% start with a distance above the convergence criterion
dist    = crit+1;
% counter
iter    = 1;

fprintf('Inner loop, running... \n');
    
    while (dist>crit&&iter<maxiter)

        % derive current consumption
        c0 = C0(cp0,r0); % cp0 is a guess of expected future consumption 

        % derive current assets
        a0 = A0(Amat,Ymat,c0,r0,w0);

        %%% update the guess for the consumption decision

        % consumption decision rule for a binding borrowing constraint
        % can be solved as a quadratic equation
        % c^(2)-((1+r)a+b)c-(yw)^(1+2/3) = 0
        % general formula and notation: ax^(2)+bx+c = 0
        % x = (-b+sqrt(b^2-4ac))/2a
        %cpbind = ((1+r0)*Amat+b)/2+sqrt(((1+r0)*Amat+b).^(2)+4*(Ymat*w0).^(1+varphi))/2;
        cpbind = (1+r0)*Amat-b+Ymat*w0;
        % % slow alternative: rootfinding gives the same result
        % cpbind  = zeros(M,N);
        % options = optimset('Display','off');
        % cpbind(:,1) = fsolve(@(c) (1+r0)*Amat(:,1)+H(c,Y(1),w0).*Y(1).*w0-c+b,cp0(:,1),options);
        % cpbind(:,2) = fsolve(@(c) (1+r0)*Amat(:,2)+H(c,Y(2),w0).*Y(2).*w0-c+b,cp0(:,2),options);

        % consumption for nonbinding borrowing constraint
        cpnon = zeros(M,N);
        % interpolation conditional on productivity realization
        % instead of extrapolation use the highest value in a0
        cpnon(:,1)  = interp1(a0(:,1),c0(:,1),Amat(:,1),'spline');
        cpnon(:,2)  = interp1(a0(:,2),c0(:,2),Amat(:,2),'spline');

        % merge the two, separate grid points that induce a binding borrowing constraint
        % for the future asset level (the first observation of the endogenous current asset grid is the 
        % threshold where the borrowing constraint starts binding, for all lower values it will also bind
        % as the future asset holdings are monotonically increasing in the current
        % asset holdings).
        cpnext(:,1) = (Amat(:,1)>a0(1,1)).*cpnon(:,1)+(Amat(:,1)<=a0(1,1)).*cpbind(:,1);
        cpnext(:,2) = (Amat(:,2)>a0(1,2)).*cpnon(:,2)+(Amat(:,2)<=a0(1,2)).*cpbind(:,2);

        % distance measure
        dist = norm((cpnext-cp0)./cp0);

        % display every 100th iteration
        if mod(iter,100) == 1
            fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dist]);
        end

        % increase counter
         iter = iter+1;

        % update the guess on the consumption function
         cp0 = cpnext;
    end
% The grid is endogenous so the asset policy function will be: 
%aprime = a0*(1+r0)+H(c0,Ymat,w0).*Ymat.*w0-c0; % Not taking into account that the borrowing constraint  is being violated
%        cpnext(:,1) = (aprime(:,1)>b).*aprime(:,1)+(aprime(:,1)<=b).*cpbind(:,1);
%        cpnext(:,2) = (aprime(:,2)>b).*aprime(:,2)+(aprime(:,2)<=b).*cpbind(:,2);
%figure(2)
%plot(a0,Amat)
% Prepare for interpolation
%
clear aprime
agrid = linspace(b,aM,M);
aprime(:,1) = interp1(a0(:,1),Amat(:,1),agrid,'pchip');
aprime(:,2) = interp1(a0(:,2),Amat(:,2),agrid,'pchip');
cpol(:,1) = interp1(a0(:,1),cp0(:,1),agrid,'pchip');
cpol(:,2) = interp1(a0(:,2),cp0(:,2),agrid,'pchip');
%a0_tream = a0(a_aux1,2);
figure(1)
plot(agrid,cpol)
title('Consumption')
figure(4)
plot(agrid',aprime)
title('Assets')
%}
%end
%M  = 100;
K = 2*M;
agrid_finer = linspace(b,aM,K);
%endog_grid = Amat;
endog_grid = a0;
Lambda0 = zeros(K,N);
for j = 1:N
    for k = 1:K 
        Lambda0(k,j) = (agrid_finer(k)-b)/(aM-b)*invdist(j);
    end
end
Lnm1_mat = Lambda0;
Ln_mat   = zeros(K, N) ;
  iter2 = 0; 
  dist = 10; 


[amesh,ymesh]=meshgrid(agrid,Y');
%
% With 3D dimensional
%a_ast_itp = @(a_ipo,b_ipo) interp2(amesh,ymesh,endog_grid',a_ipo,b_ipo,'spline'); % whi not linear, cuz NaN will always be ataken as a minimum

% With 2D dimensional
%a_ast_itp = @(a_ipo,pos) interp1(agrid,endog_grid(:,pos),a_ipo,'linear','extrap');
a_ast_itp = @(a_ipo,pos) interp1(agrid,endog_grid(:,pos),a_ipo,'pchip'); 
% Make the distribution contonous
% Make it continous % This is the mistake
%l0 = @(a,lambdarg) min(max(interp1(agrid_finer,lambdarg,a,'linear','extrap'),0),1); % not linear cuz it doesnt extrapolater
%mesh(amesh,ymesh,endog_grid') 
X0 = 0;
options = optimset('Display','off');
    while iter2<250 && dist>1e-4 
 
       iter2 = iter2+1 

 
      for i_y = 1:N %# next period 
           for i_a =1:K %# next period 
           %for i_a =K:1 %# next period 
              % a_v is the value
                    % i_a is the index
              a_v = agrid_finer(i_a);
              vec_temp = zeros(N,1); 

                 for i_y0 = 1:N   % # last period y 
                     y0_v = Y(i_y0);
                     %3 D
                     %aval  = min(max(a_ast_itp(a_v, y0_v), agrid(1)), agrid(end));% # today's assets (endogenous grid) 
                     %2 D
                     aval  = min(max(a_ast_itp(a_v, i_y0), agrid(1)), agrid(end));% # today's assets (endogenous grid) 
                     % Solver
                     %aval = 
                     %aval_aux=fsolve(@(a) invaprime(a,Y,r0,a0,cp0,agrid_finer(i_a),i_y0,w0),X0,options);
                     %aval = min(max(aval_aux,agrid(1)),agrid(end));
                     
                     %ind_r = min(max(searchsortedfirst(ay.agrid_finer, aval), 2), K) 
                     %index_aux = nearestpoint(aval,agrid_finer,'nearest'); % nice trick
                     index_aux = find(agrid_finer>=aval,1,'first');
                     if isempty(index_aux)
                          index_aux=K;
                     end
                     %ind_r = index_aux;
                     %ind_r = min(index_aux,K-1);
                     ind_r = max(index_aux,2);
 %                   ?val = ?nm1_mat[ind_r-1,i_y0] + (?nm1_mat[ind_r, i_y0]- ?nm1_mat[ind_r-1, i_y0]) / (agrid_finer[ind_r] - agrid_finer[ind_r-1]) * (aval - agrid_finer[ind_r-1]) 
                     Lval = Lnm1_mat(ind_r-1,i_y0) + (Lnm1_mat(ind_r, i_y0)- Lnm1_mat(ind_r-1, i_y0)) / (agrid_finer(ind_r) - agrid_finer(ind_r-1)) * (aval - agrid_finer(ind_r-1)); 
                     %Lval = Lnm1_mat(ind_r,i_y0) + (Lnm1_mat(ind_r+1, i_y0)- Lnm1_mat(ind_r, i_y0)) / (agrid_finer(ind_r+1) - agrid_finer(ind_r)) * (aval - agrid_finer(ind_r)); 
                     %Lval = Lnm1_mat(ind_r,i_y0) + l0(aval,Lnm1_mat(:,i_y));
                     vec_temp(i_y0) = Lval; 
                 end 
 
 
                 %Ln_mat(i_a, i_y) = dot(pi(:, i_y),vec_temp); 
                  Ln_mat(i_a, i_y) = dot(pi(:, i_y),vec_temp); 
            end 
      end 
         dist = max(max(abs(Ln_mat - Lnm1_mat))); 
         Lnm1_mat=Ln_mat;
         %copy!(?nm1_mat, ?n_mat) 

         %# if iter%200 == 0.0 
         %#     @printf("Iteration # %d on distribution with distance %.3g\n", iter, dist) 
         %# end 
    end 
    
    invariant = Lnm1_mat;
    %copy!(?_invariant, ?nm1_mat) 
   % Void 
 
pr = sum(invariant,2)';
pr0 = sum(Lambda0,2)';
figure(10)
plot(agrid_finer,pr);
figure(11)
plot(agrid_finer,pr0);
%
grid_even_finer = linspace(b,aM,300*K);
new_pr = interp1(agrid_finer,pr,grid_even_finer,'pchip');

normalization = 1-new_pr(end);
newpr= new_pr+normalization;
pdfpr = diff(newpr);
figure(12)
plot(grid_even_finer,[0,pdfpr]);
%}

%% transition matrix
%{
% Piecewise linear interpolation of the invariant distribution
%b = -2; aM = 20; %M  = 200; K=M;
K = 100;%2000;
agrid2 = linspace(b,aM,K);
Lambda0 = zeros(K,N);
for j = 1:N
    for k = 1:K 
        Lambda0(k,j) = (agrid2(k)-b)/(aM-b)*invdist(j);
    end
end
% Make it continous % This is the mistake
l0 = @(a,lambdarg) interp1(agrid2,lambdarg,a,'pchip'); % not linear cuz it doesnt extrapolater
% update distribution for every pair
crit1 = 0.001;%10^(-10);
dist1 = 1;
lambda = Lambda0;
l1 = zeros(K,N);
X0 = 0;
iter2 = 0;
options = optimset('Display','off');
%
while dist1>crit1 && iter2<130%dist1>crit1 %dist>crit&&iter<maxiter
    iter2 = iter2+1    
        for k = 1:K
            for j = 1:N
                for i = 1:N
                    
                    %fsolve(@(x)sin(x.*x),x0);
                    asol_aux=fsolve(@(a) invaprime(a,Y,r0,a0,cp0,agrid2(k),i,w0),X0,options);
                    asol = min(max(asol_aux,agrid(1)),agrid(end));
                    l1(k,j) = l1(k,j)+pi(j,i)*l0(asol,lambda(:,i));
                    %l1(k,j) = l1(k,j)+pi(j,i)*Linterp(k,i,lambda,asol,agrid2);
                    dert = 1;
                end
            end
        end

    l2 = l1;
    dist1 = max(max(abs(lambda-l2)));
    lambda = l1;
    l1 = zeros(K,N);
end
pr = sum(lambda,2)';
pr0 = sum(Lambda0,2)';
figure(10)
plot(agrid2,pr);
figure(11)
plot(agrid2,pr0);

normalization = 1-pr(end);
newpr= pr+normalization;
pdfpr = diff(newpr);
figure(12)
plot(agrid2,[0,pdfpr]);

%}


%asol=fsolve(@(a) invaprime(a,Y,r0,a0,cp0,agrid2(50),1),0.3);

%{





%% Discretization Transition matrix
% Creating tdecis
%[agrid',aprime];
tdecis1(:,1) =nearestpoint(aprime(:,1),agrid','nearest');
tdecis1(:,2) =nearestpoint(aprime(:,2),agrid','nearest'); 
%
gmat=zeros(M,M,N);
   
   for j = 1:N
      
      for k = 1:M
         
         gmat(k,tdecis1(k,j),j) = 1;
         
      end
      
      trans((j-1)*M+1:j*M,:) = kron(pi(j,:),gmat(:,:,j));
         
   end
   
   trans=trans';
   
   probst = (1/(N*M))*ones(N*M,1);
   test=1;
   while test > 10^(-5);
      probst1 = trans*probst;
      test = max(abs(probst1-probst));
      probst = probst1;
   end;
   %
   %   vectorize the decision rule to be conformable with probst
   %   calculate new aggregate capital stock  meanK
   %
   
   kk=aprime(:);
   meank=probst'*kk
   
figure(2)
pp = reshape(probst,[length(probst)/N N]);
pr = sum(pp,2)';
plot(agrid',pr);
%
p = polyfit(agrid,pr,11);
x1 = linspace(b+0.01,aM,500);
y1 = polyval(p,x1);
%
figure(4)
plot(x1,y1)
%}