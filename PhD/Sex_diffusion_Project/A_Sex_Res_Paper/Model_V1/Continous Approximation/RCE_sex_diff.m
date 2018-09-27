clear all
% Load Parameters
[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(1);
%{
% Solver Options
options_2 = optimset('Display','iter');
options_1 = optimset('Display','iter','MaxFunEvals',1e2,'MaxIter',1e2,'TolFun',1e-10,'TolX',1e-5);

% Initial Values
rate_init = 0.02;
price_init = 0.6;
rate_price = [rate_init,price_init];
type      = 1;
gender    = 1;
% Solution interest rate
%[price,sex]  = fsolve(@(price) sex_agg(rate_init,price),price_init,options_1) 
%[rate,asset] = fzero(@(zero_rate) excess(zero_rate,type),rate_init,options_2)
%[rate,asset] = fsolve(@(zero_rate) a_excess(zero_rate,1,method),rate_init)
%}
%% Choose a type
type = 1;
%% solve 
%the limits and the initial mesh
ax=[];
ax(1).val=linspace(0,1,20);
ax(2).val=linspace(eps,100,20);

% number of iteration
Niteration=4;% ori7take care, the large values can easily lead to memory problem

%% Bisection solver
%{
bound_fuction_name = 'sol_aux';
outputsol  = db_mdbm(ax,bound_fuction_name,Niteration);
%% Retrieve solution
d_solution = outputsol.posinterp;
rate  = d_solution(1);
price = d_solution(2);
rate_price = [rate,price];
%}

%% Fsolve solver
%{
options_1 = optimset('Display','iter');
[rate_price] = fsolve(@(zero) clearing(zero),[0.001,0.05],options_1);
%}
% Retrive distributions
clear all
rate  = 0.0145; %  0.015/-8; 0.017/1.9;
price = 0.4;   %
rate_price = [rate,price];
[excess_total,gk]          = clearing(rate_price);

[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,~,alpha]=parameters(1);
figure(9)
plot(ag',sum(gk,2));
title('Asset distribution');
xlabel('Asset holdings');
ylabel('% OF AGENTS');
%%
rate  = 0.08; %  0.015/-8; 0.017/1.9;
price = 0.32;
[aopt_m_edu,copt_m_edu,xopt_m_edu,v_m_edu] = partial_gs(rate,price,1,1);
[aopt_f_edu,copt_f_edu,xopt_f_edu,v_f_edu] = partial_gs(rate,price,1,2);
[aopt_m_nedu,copt_m_nedu,xopt_m_nedu,v_m_nedu] = partial_gs(rate,price,2,1);
[aopt_f_nedu,copt_f_nedu,xopt_f_nedu,v_f_nedu] = partial_gs(rate,price,2,2);
%}
%
%% save some results
%{
Try1 with 100, degen, eps,100, 20. Iterval -5 8. distribution swiked to the
left, No market lcearing in the result. iterations over ergo, 5000*4;
Try2 do it with the fsolve rate:0.0047 price:0.6778. excess: -7.1 2.9


%}

%% Figures
%clear all
[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(1);
%[aopt_m,copt_m,xopt_m,v_m] = partial_male(0.0035,0.056,1);
%[aopt_m,copt_m,xopt_m,v_m] = partial_male_interp(0.0035,0.056,1);
%[aopt_m_edu,copt_m_edu,xopt_m_edu,v_m_edu] = partial_gs(0.02,0.6,1,1); %-19.7
%[aopt_m_edu,copt_m_edu,xopt_m_edu,v_m_edu] = partial_gs(0.015,0.8,1,1); %-19.7
%[aopt_m_edu,copt_m_edu,xopt_m_edu,v_m_edu] = partial_gs(0.08,0.32,1,1); %-19.7 %new 0.185
figure(1)
plot(a',v_m_edu,'-o')
ylabel('Value Function')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Value Function Male-Educated') 
 
figure(2)
plot(a',copt_m_edu,'-o')
ylabel('Consumption')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Consumption Male-Educated')

figure(3)
plot(a',aopt_m_edu,'-o')
ylabel('Next-Period asset holdings')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Assets Male-Educated')

figure(4)
plot(a',xopt_m_edu,'-o')
ylabel('Sex consumption')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Sex consumption Male-Educated')
beep()
%%
%[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(1);
figure(5)
plot(a',v_f_edu,'-o')
ylabel('Value Function')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Value Function Female-Educated') 
 
figure(6)
plot(a',copt_f_edu,'-o')
ylabel('Consumption')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Consumption Female-Educated')

figure(7)
plot(a',aopt_f_edu,'-o')
ylabel('Next-Period asset holdings')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Assets Female-Educated')

figure(8)
plot(a',xopt_f_edu,'-o')
ylabel('Sex production')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Sex Production Female-Educated')
%}
%%
[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,~,alpha]=parameters(1);
figure(9)
plot(ag',sum(gk,2));
title('Asset distribution');
xlabel('Asset holdings');
ylabel('% OF AGENTS');
%}
%% Retrieving income distrution 

[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st,~,alpha]=parameters(1);

inc_m_edu      = zeros(size(ag,2),n_st);
noc_inc_m_edu  = zeros(size(ag,2),n_st);
inc_f_edu      = zeros(size(ag,2),n_st);
noc_inc_f_edu  = zeros(size(ag,2),n_st);
inc_m_nedu     = zeros(size(ag,2),n_st);
noc_inc_m_nedu = zeros(size(ag,2),n_st);
inc_f_nedu     = zeros(size(ag,2),n_st);
noc_inc_f_nedu = zeros(size(ag,2),n_st);
for i =1:n_st 
    % Educated
    % Buyers % Addinterpolations to increase the number of points
    inc_m_edu(:,i)      = endow(i,1)+(1+rate).*(ag');  %loop for the male_edus
    noc_inc_m_edu(:,i)  = interp1(a',copt_m_edu(:,i),ag','spline') + (price.*interp1(a',xopt_m_edu(:,i),ag','spline'))+interp1(a',aopt_m_edu(:,i),ag','spline')-(1+rate).*(ag');
                         
    % Sellers 
    lab_f_edu           = (endow(i,1)/(alpha*price))^(1/(alpha-1)); sex_f_edu = lab_f_edu^alpha;
    inc_f_edu(:,i)      = endow(i,1)*(1-lab_f_edu)+(1+rate).*(ag')+price*sex_f_edu;
    noc_inc_f_edu(:,i)  = interp1(a',copt_f_edu(:,i),ag','spline')+interp1(a',aopt_f_edu(:,i),ag','spline')-(1+rate).*(ag');
    
    % Less Educated
    % Buyers
    inc_m_nedu(:,i)     = endow(i,2)+(1+rate).*(ag');  %loop for male_nedus
    noc_inc_m_nedu(:,i) = interp1(a',copt_m_nedu(:,i),ag','spline') + (price.*interp1(a',xopt_m_nedu(:,i),ag','spline'))+interp1(a',aopt_m_nedu(:,i),ag','spline')-(1+rate).*(ag');
    % Sellers
    lab_f_nedu          = (endow(i,2)/(alpha*price))^(1/(alpha-1)); sex_f_nedu = lab_f_nedu^alpha;
    inc_f_nedu(:,i)     = endow(i,2)*(1+lab_f_nedu)+(1+rate).*(ag')+price*sex_f_nedu;
    noc_inc_f_nedu(:,i) = interp1(a',copt_f_nedu(:,i),ag','spline')+interp1(a',aopt_f_nedu(:,i),ag','spline')-(1+rate).*(ag');
end

%%
% For graph 13
figure(13)
plot(ag',noc_inc_m_edu(:,1),ag',noc_inc_m_edu(:,2),ag',noc_inc_f_edu(:,1),ag',noc_inc_f_edu(:,2));
title('No capital income (NCI)');
xlabel('Asset holdings');
ylabel('NCI');
legend('Buyers-G','Buyers-B','Sellers-G','Sellers-B')
%%
%dist = [d_m_edu,d_f_edu,d_m_nedu,d_f_nedu]; %16 column
% create the income matrix
inc = [inc_m_edu,inc_f_edu,inc_m_nedu,inc_f_nedu]; %16 colums
aux_inc = inc(:);
min_inc = min(aux_inc); % Should not be negative
max_inc = max(aux_inc);
income_grid = linspace(min_inc,max_inc,size(ag,2));
income_dist = zeros(1,size(ag,2));
for i = 1:n_st*n_st
    for j = 1:size(ag,2)
    [~,index_inc] = min(abs(income_grid-inc(j,i))); % Finds the closesnt point in the grid to inc(j,i)
    income_dist(1,index_inc) = income_dist(1,index_inc)+gk(j,i); % gk 16 columns
    end
end
figure(10)
plot(income_grid,income_dist);
title('Income distribution');
xlabel('Income == Y+(1+r)a');
ylabel('% OF AGENTS');
%% Combined graphs
[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(1);
figure(11)
plot(a',aopt_m_edu(:,1),a',aopt_m_edu(:,2),a',aopt_f_edu(:,1),a',aopt_f_edu(:,2));
title('Asset policy function');
xlabel('Asset holdings');
ylabel('Next period asset holdings');
legend('Buyers-G','Buyers-B','Sellers-G','Sellers-B')

figure(12)
plot(a',copt_m_edu(:,1),a',copt_m_edu(:,2),a',copt_f_edu(:,1),a',copt_f_edu(:,2));
title('Consumption policy function');
xlabel('Asset holdings');
ylabel('Consumption');
legend('Buyers-G','Buyers-B','Sellers-G','Sellers-B')
%%
%
for k=1:8
  saveas(figure(k),sprintf('FIG%d.png',k)); % will create FIG1, FIG2,...
end
%}