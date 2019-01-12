function [aopt,copt,xopt,v]=partial_gs(r,price,type,gender)
%[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(p)
[sigma,endow_all,~,a,na,beta,pp,eh_g,el_g,pp1,amin1,amax1,~,~,~,~,n_st,ngk,alpha,nit,tol,epsilon,tol1,gamma,phy]=parameters(1);
if type == 1
    endow = endow_all(:,type);
else
    endow = endow_all(:,type);
end

             % /* agents consume their income */
copt = zeros(na,n_st);          % /* optimal consumption */
%aopt = zeros(na,n_st);          % /* optimal next-period assets */
xopt = zeros(na,n_st);          % /* optimal next-period assets */
v0 = zeros(na,n_st);
v1 = zeros(na,n_st);
aopt = zeros(na,n_st);
dif = 1;
its = 0;
while dif > 10e-10
    %while dif>tol && its < maxits
        for i = 1:na
            for j = 1:n_st
            a0 = a(1,i);
            a1 = fminbnd(@(a_prime)bellman_int(a_prime,j,v0,r,price,type,a0,gender),amin1,amax1);
            v1(i,j) = -bellman_int(a1,j,v0,r,price,type,a0,gender);
            aopt(i,j) = a1;
            end
        end
    dif = norm(v1-v0);
    v0 = v1;
    its = its+1;
end
v = v0;
% retrive policy functions
if gender == 1
		copt(:,1) = (endow(1) + (1+r).*a'-aopt(:,1))./(1+price.^(1-(1/sigma)));
        xopt(:,1) = copt(:,1).*price.^(1/-sigma);
		copt(:,2) = (endow(2) + (1+r).*a'-aopt(:,2))./(1+price.^(1-(1/sigma)));
        xopt(:,2) = copt(:,2).*price.^(1/-sigma);
		copt(:,3) = (endow(3) + (1+r).*a'-aopt(:,3))./(1+price.^(1-(1/sigma)));
        xopt(:,3) = copt(:,3).*price.^(1/-sigma);
		copt(:,4) = (endow(4) + (1+r).*a'-aopt(:,4))./(1+price.^(1-(1/sigma)));
        xopt(:,4) = copt(:,4).*price.^(1/-sigma);
else
      lopt(:,1) = repmat((endow(1)/(alpha*price))^(1/(alpha-1)),na,1);
      xopt(:,1) = lopt(:,1).^alpha;
      copt(:,1) = price*xopt(:,1)+(endow(1)*(1-lopt(:,1)))+(1+r).*a'-aopt(:,1);

	  lopt(:,2) = repmat((endow(2)/(alpha*price))^(1/(alpha-1)),na,1);
      xopt(:,2) = lopt(:,2).^alpha;
      copt(:,2) = price*xopt(:,2)+(endow(2)*(1-lopt(:,2)))+(1+r).*a'-aopt(:,2);
	  
	  lopt(:,3) = repmat((endow(3)/(alpha*price))^(1/(alpha-1)),na,1);
      xopt(:,3) = lopt(:,3).^alpha;
      copt(:,3) = price*xopt(:,3)+(endow(3)*(1-lopt(:,3)))+(1+r).*a'-aopt(:,3);
	  
	  lopt(:,4) = repmat((endow(4)/(alpha*price))^(1/(alpha-1)),na,1);
      xopt(:,4) = lopt(:,4).^alpha;
      copt(:,4) = price*xopt(:,4)+(endow(4)*(1-lopt(:,4)))+(1+r).*a'-aopt(:,4);
end


end