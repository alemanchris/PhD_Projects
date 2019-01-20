function [aopt,copt,xopt,v]=partial_gs(r,price,type,gender)
%[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(p)
[sigma,endow_all,~,a,na,beta,pp,eh_g,el_g,pp1,amin1,amax1,~,~,~,~,n_st,ngk,alpha,nit,tol,epsilon,tol1,gamma_all,phy]=parameters(1);
 if type==1 
     endow = endow_all(:,1);
     %gamma = gamma_all(:,1);
 else
     endow = endow_all(:,2);
     %gamma = gamma_all(:,2);
 end
             
copt = zeros(na,n_st);          % /* optimal consumption */
xopt = zeros(na,n_st);          % /* optimal sex-consumption */
v0 = zeros(na,n_st);
v1 = zeros(na,n_st);
aopt = zeros(na,n_st);
dif = 1;
its = 0;
%while dif > 10e-10
maxits = 1000;
 while dif>tol && its < maxits
        for i = 1:na
            for j = 1:n_st
            a0 = a(1,i);
            a1 = fminbnd(@(a_prime)bellman_int(a_prime,j,v0,r,price,type,a0,gender),amin1-10e-6,amax1+10e-6); % I add this cuz the maximization does not include the bounds
            %fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
            %a1 = fminbnd(@(a_prime)bellman_int(a_prime,j,v0,r,price,type,a0,gender),amin1,amax1);
            %a1 = fmincon(@(a_prime)bellman_int(a_prime,j,v0,r,price,type,a0,gender),a0,[],[],[],[],amin1,amax1);

            v1(i,j) = -bellman_int(a1,j,v0,r,price,type,a0,gender);
            aopt(i,j) = a1;
            end
        end
    %dif = norm(v1-v0);
    dif =max(max(abs((v1-v0)./v1)));
    v0  = v1;
    its = its+1;
    disp(its)
end
v = v0;
% retrive policy functions
if gender ==1
    for i =1:n_st
            copt(:,i) = (endow(i) + (1+r).*a'-aopt(:,i))./(1+price.^(1-(1/sigma)));
            xopt(:,i) = copt(:,i).*price.^(1/-sigma);
    end
else
        lopt = zeros(1,n_st);
        for i =1:n_st
            lopt(i)   = (endow(i)/(alpha*price))^(1/(alpha-1));
            xopt(:,i) = repmat(lopt(i)^alpha,na,1);
            copt(:,i) = price*xopt(:,i)+(endow(i)*(1-lopt(i)))+(1+r)*a'-aopt(:,i);
	  
        end
end
%{
		copt(:,1) = (endow(1) + (1+r).*a'-aopt(:,1))./(1+price.^(1-(1/sigma)));
        xopt(:,1) = copt(:,1).*price.^(1/-sigma);
		copt(:,2) = (endow(2) + (1+r).*a'-aopt(:,2))./(1+price.^(1-(1/sigma)));
        xopt(:,2) = copt(:,2).*price.^(1/-sigma);
		copt(:,3) = (endow(3) + (1+r).*a'-aopt(:,3))./(1+price.^(1-(1/sigma)));
        xopt(:,3) = copt(:,3).*price.^(1/-sigma);
		copt(:,4) = (endow(4) + (1+r).*a'-aopt(:,4))./(1+price.^(1-(1/sigma)));
        xopt(:,4) = copt(:,4).*price.^(1/-sigma);
 %}
end