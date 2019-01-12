function val=bellman_int(k,j,v0,r,price,type,k0,gender)
% This program gets the value function for a neoclassical growth model with
% no uncertainty and CRRA utility
[sigma,endow_all,~,a,na,beta,pp,~,~,~,~,~,~,~,~,~,n_st,ngk,alpha,nit,tol,epsilon,tol1,gamma_all,phy]=parameters(1);

%global v0 beta delta alpha kmat k0 s kgrid y P
if type == 1
    endow = endow_all(:,type);
    gamma = gamma_all(:,type);
else
    endow = endow_all(:,type);
    gamma = gamma_all(:,type);
end

klo=max(sum(k>a),1); % identify the gridpoint that falls just below . .
 % . . the choice for k
 khi=klo+1;

% do the interpolation
%gg = v0(klo,:) + (k-a(klo))*(v0(khi,:) - v0(klo,:))/(a(khi) - a(klo));
%gg = v0(klo,) + (k-a(klo))*(v0(khi,j) - v0(klo,j))/(a(khi) - a(klo));
gg = interp1(a,v0,k,'spline');
%interp1(a,v,kp,'spline'
 %c = k0^alpha - k + (1-delta)*k0; % consumption
% c = k0*(1+r) + endow(j) - k; % consumption
if gender==1
     c = (endow(j) + (1+r)*k0-k)/(1+price^(1-(1/sigma)));
     x = c*price^(1/-sigma);
else
     l = (endow(j)/(alpha*price))^(1/(alpha-1));
     x = l^alpha;
     c = price*x+(endow(j)*(1-l))+(1+r)*k0-k;
end
if c<=0 || x<=0
val = -9999999 - 999*abs(c);
else
%val = (1/(1-s))*(c^(1-s)) + beta*((gg(1,1)*P(j,1))+(gg(1,2)*P(j,2)));
val = u(c,x,gender) + beta*gamma(j,1)*((gg(1,1)*pp(j,1))+(gg(1,2)*pp(j,2))+(gg(1,3)*pp(j,3))+(gg(1,4)*pp(j,4)));
 end
 val = -val; 