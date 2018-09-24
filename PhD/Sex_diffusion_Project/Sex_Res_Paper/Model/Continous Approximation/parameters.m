function [sigma,endow,neg,a,na,beta,pp,eh_g,el_g,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,ngk,alpha,nit,tol,epsilon,tol1,gamma,phy]=parameters(p)
%[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,nkg,alpha]=parameters(1);
if p==1
%/* parameter values */
beta    =0.99322;
alpha   =0.1;
sigma	=1.5;
neg		=-1e10;

%/* fertility values if Onset*/
gamma    =[0.98,0.97; %hg hb lg lb
           0.98,0.97;
           0.90,0.89;
           0.90,0.89]; % Survival rate = gamma \ Mortality rate (1-gamma)
phy    =[0.7,0.7,0.4,0.4; %Fertility rate Educated
         0.8,0.8,0.5,0.5];% 
     
     %/* fertility values if Pre-Epidemic*/
gamma    =[0.98,0.97; %hg hb lg lb
           0.98,0.97;
           0.98,0.97;
           0.98,0.97]; % Survival rate = gamma \ Mortality rate (1-gamma)
phy    =[0.7,0.7,0.7,0.7; %Fertility rate Educated
         0.8,0.8,0.8,0.8];% 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %{
		      %/* fertility values if Pre-Epidemic*/
gamma    =[1; %hg hb lg lb
           1;
           1;
           1]; % Survival rate = gamma \ Mortality rate (1-gamma)
phy    =[0,0,0,0; %Fertility rate Educated
         0,0,0,0];% 
     %}


%/* endowment Onset*/
eh_g      = 1;			% Healthy_goodtimes
el_g	  = 0.7;		% Infected_goodtimes
eh_b      = 0.4;        % Healthy_badtimes
el_b      = 0.1;        % Infected_badtimes
%
%/* endowment if preepidemic */
eh_g	  = 1;			% Healthy_goodtimes
el_g	  = 1;			% Infected_goodtimes
eh_b      = 0.4;          % Healthy_badtimes
el_b      = 0.4;          % Infected_badtimes
% Education heterogeneity
educ = 0.2;
%}
endow	= [eh_g,eh_g-educ;
           eh_b,eh_b-educ;
           el_g,el_g-educ;
           el_b,el_b-educ];  %second column is non edu
%pp		= [0.925 0.075;0.5 0.5];
prob_tauchen_h   = [0.70 0.30;
                    0.40 0.60]; % prob(i,j) = probability (s(t+1)=sj | s(t) = si)
prob_infection   = [0.50 0.50;
                    0.50 0.50];
%prob_inf   = [ .8 .2; 0.4 0.6]; % prob(i,j) = probability (s(t+1)=sj | s(t) = si)

pp = kron(prob_infection,prob_tauchen_h);% makes hg hb lg lb
%/* ergodic distribution */
pp1		= equivec1(pp);
% /* asset grid */
amin1	= -5;   % -2 works well with-5  15   9         
amax1	= 28;    %  4 works well with 8  28   15
na		= 7;  %/*na=101;*/ na30 300 funca 70 ultimo % 20
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


%{
function [sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,ngk,alpha]=parameters(p)
%[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,nkg,alpha]=parameters(1);
if p==1
%/* parameter values */
beta    =0.99322;
alpha   =0.1;
sigma	=1.5;
neg		=-1e10;

%/* endowment if not preepidemic */
eh		= 1;			% Healthy_good
el		= 0.7;			% Healthy_bad
ej      = 0.4;
eq      = 0.1;
%/* endowment if preepidemic */
eh		= 1;			% Healthy_good
el		= 0.4;			% Healthy_bad
ej      = 1;
eq      = 0.4;
endow	= [eh,eh; %Keep in mind the order is wrong, the order in continous is correct
		   el,el;
		   ej,ej;
		   eq,eq];  %second column is non edu
%pp		= [0.925 0.075;0.5 0.5];
prob_tau   = [ .65 .35; .5 .5]; % prob(i,j) = probability (s(t+1)=sj | s(t) = si)
prob_inf   = [ 0.5 0.5; 0.5 0.5];
%prob_inf   = [ .8 .2; 0.4 0.6]; % prob(i,j) = probability (s(t+1)=sj | s(t) = si)

pp = kron(prob_inf,prob_tau);% makes hg hu ie iu
%/* ergodic distribution */
pp1		= equivec1(pp);
% /* asset grid */
amin1	= -2;   % -2 works well with-5  15   9         
amax1	= 8;    %  4 works well with 8  28   15
na		= 70;  %/*na=101;*/ na30 300 funca
astep	= (amax1-amin1)/(na-1);
a		= linspace(amin1,amax1,na);
%/* asset grid for distribution */
%nk		= 3*na;   % 3 times larger   3*na        
nk		= 2*na;
ngk		= 500*5;         % /* 5000*5 number of iterations over distribution 5000*/
agstep	= (amax1-amin1)/(nk-1);
ag		= linspace(amin1,amax1,nk);
n_st    = 4;
end


end
%}