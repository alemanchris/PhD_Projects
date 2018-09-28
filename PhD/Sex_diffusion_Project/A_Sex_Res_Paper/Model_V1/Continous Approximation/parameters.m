function [sigma,endow,neg,a,na,beta,pp,eh_g,el_g,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,ngk,alpha,nit,tol,epsilon,tol1,gamma,phy]=parameters(p)
%[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,nkg,alpha]=parameters(1);
if p==1
%/* parameter values */
beta    =0.99322;
beta    = 0.9;
alpha   =0.1;
sigma	=1.5;
neg		=-1e10;

%/* fertility values if Onset*/
%{
gamma    =[0.99,0.985; %hg hb lg lb worked with 98 all
           0.99,0.985;
           0.98,0.975;
           0.98,0.975]; % Survival rate = gamma \ Mortality rate (1-gamma)
       %
phy    =[0.7,0.7,0.4,0.4; %Fertility rate Educated
         0.8,0.8,0.5,0.5];% 
%
%/* endowment Onset*/ % makes hg hb lg lb
eh_g      = 0.8;	    % Healthy_goodtimes  8
el_g	  = 0.4;		% Infected_goodtimes 4
eh_b      = 0.75;        % Healthy_badtimes  7
el_b      = 0.35;        % Infected_badtimes 3 
%
%Transition Onset
prob_tauchen_h   = [0.95 0.05;
                    0.05 0.95]; % prob(i,j) = probability (s(t+1)=sj | s(t) = si)
prob_infection   = [0.50 0.50;
                    0.00 1.00];
 %}    
%/* fertility values if Pre-Epidemic*/
%
gamma    =[0.99,0.98; %hg hb lg lb
           0.99,0.98;
           0.99,0.98;
           0.99,0.98]; % Survival rate = gamma \ Mortality rate (1-gamma)
       %
phy    =[0.02,0.02,0.02,0.02; %Fertility rate Educated
         0.03,0.03,0.03,0.03];% 
%phy    =[0,0,0,0; %Fertility rate Educated
%         0,0,0,0];% 
     shock = 0.8; 
     % 0.6;
%/* endowment if preepidemic */
eh_g	  = 1;		%0.8	% Healthy_goodtimes
el_g	  = 1;		%0.8	% Infected_goodtimes
eh_b      = eh_g*(1-shock);          % Healthy_badtimes 0.4 0.7
el_b      = el_g*(1-shock);          % Infected_badtimes

% Transition pree-epidemic
%prob_tauchen_h   = [0.70 0.30;
%                    0.40 0.60];
prob_tauchen_h   = [0.55 0.45;
                    0.45 0.55];% prob(i,j) = probability (s(t+1)=sj | s(t) = si)
prob_infection   = [0.50 0.50;
                    0.50 0.50];
%}
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %
		      %/* fertility values if Pre-Epidemic*/
     %{
gamma    =[1,1; %hg hb lg lb
           1,1;
           1,1;
           1,1]; % Survival rate = gamma \ Mortality rate (1-gamma)
       %
phy    =[0,0,0,0; %Fertility rate Educated
         0,0,0,0];% 
     %}
%}



% Education heterogeneity
educ = 0.05;

%}
endow	= [eh_g,eh_g*(1-educ);
           eh_b,eh_b*(1-educ);
           el_g,el_g*(1-educ);
           el_b,el_b*(1-educ)];  %second column is non edu



pp = kron(prob_infection,prob_tauchen_h);% makes hg hb lg lb
%/* ergodic distribution */
pp1		= equivec1(pp);
% /* asset grid */
%+10e-6
amin1	= -10;   % -2 works well with-5  15   9         
amax1	= 10;    %  4 works well with 8  28   15
na		= 12;  %/*na=101;*/ na30 300 funca 70 ultimo % 20
astep	= (amax1-amin1)/(na-1);
a		= linspace(amin1,amax1,na);
%/* asset grid for distribution */
%nk		= 3*na;   % 3 times larger   3*na        
nk		= 100;
ngk		= 500*5;         % /* 5000*5 number of iterations over distribution 5000*/
agstep	= (amax1-amin1)/(nk-1);
ag		= linspace(amin1,amax1,nk);
n_st    = 4;

% % Computational parameters 
tolw 	=1e-6;       %  /* stopping criterion for value function */
tol  	=0.001;      %  /* stopping criterion for absolute divergence in capital stock */
tol1 	=1e-7;       %  /* stopping criterion for golden section search */
tolg	=1e-10;      %  /* distribution function */
delr	=0.01;       %  /* if excess demand is 1, the interest rate is decreased by zeta */
epsilon =0.05;
nit		=250;          % /* number of iterations over value function 250*/

nq 		=100;          % /* number of iterations over r 100*/                   
psi		=0.5;

end


end


