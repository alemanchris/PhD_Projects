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
rate  = 0.0035; % 0.006
price = 0.056;   %rate: 0.0055 price: 0.0566
rate_price = [rate,price];
[excess_total,gk]          = clearing(rate_price);
%rate  = rate_price(1);
%price = rate_price(2);
type =1;
[aopt_m,copt_m,xopt_m,v_m] = partial_gs(rate,price,type,1);
[aopt_f,copt_f,xopt_f,v_f] = partial_gs(rate,price,type,2);
%
%% save some results
%{
Try1 with 100, degen, eps,100, 20. Iterval -5 8. distribution swiked to the
left, No market lcearing in the result. iterations over ergo, 5000*4;
Try2 do it with the fsolve rate:0.0047 price:0.6778. excess: -7.1 2.9


%}

%% Figures
clear all
[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(1);
%[aopt_m,copt_m,xopt_m,v_m] = partial_male(0.0035,0.056,1);
%[aopt_m,copt_m,xopt_m,v_m] = partial_male_interp(0.0035,0.056,1);
[aopt_m,copt_m,xopt_m,v_m] = partial_gs(0.0035,0.056,1,1);
figure(1)
plot(a',v_m,'-o')
ylabel('Value Function')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Value Function Male-Educated') 
 
figure(2)
plot(a',copt_m,'-o')
ylabel('Consumption')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Consumption Male-Educated')

figure(3)
plot(a',aopt_m,'-o')
ylabel('Next-Period asset holdings')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Assets Male-Educated')

figure(4)
plot(a',xopt_m,'-o')
ylabel('Sex consumption')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Sex consumption Male-Educated')

%%
[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(1);
figure(5)
plot(a',v_f,'-o')
ylabel('Value Function')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Value Function Female-Educated') 
 
figure(6)
plot(a',copt_f,'-o')
ylabel('Consumption')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Consumption Female-Educated')

figure(7)
plot(a',aopt_f,'-o')
ylabel('Next-Period asset holdings')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Assets Female-Educated')

figure(8)
plot(a',xopt_f,'-o')
ylabel('Sex production')
xlabel('Current asset holdings')
legend('Not-Infected-G','Not-Infected_B','Infected-G','Infected_B')
title('Sex Production Female-Educated')
%}
%
figure(9)
plot(ag',sum(gk,2));
title('Asset distribution');
xlabel('Asset holdings');
ylabel('% OF AGENTS');
%}
%% Retrieving income distrution 
%% Retrieving income distrution 
inc = zeros(size(ag,n_st),n_st*2);
for i =1:n_st*2
    inc(:,i) = endow(i)+(1+rate)*(ag');
end
aux_inc = inc(:);
min_inc = min(aux_inc);
max_inc = max(aux_inc);
income_grid = linspace(min_inc,max_inc,size(ag,2));
income_dist = zeros(1,size(ag,2));
for i = 1:n_st*2
    for j = 1:size(ag,2)
    [~,index_inc] = min(abs(income_grid-inc(j,i)));
    income_dist(1,index_inc) = income_dist(1,index_inc)+gk(j,i); % gk 16 columns
    end
end
figure(10)
plot(income_grid,income_dist);
title('Income distribution');
xlabel('Income == Y+(1+r)a');
ylabel('% OF AGENTS');
%%
%
for k=1:10
  saveas(figure(k),sprintf('FIG%d.png',k)); % will create FIG1, FIG2,...
end
%}