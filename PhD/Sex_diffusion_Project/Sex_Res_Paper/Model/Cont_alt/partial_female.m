function [aopt,copt,xopt,v]=partial_female(r,price,type)
%[sigma,endow,neg,a,na,beta,pp,eh,el,pp1,amin1,amax1,astep,nk,agstep,ag,n_st]=parameters(p)
[~,endow_all,~,a,na,beta,pp,~,~,~,~,~,~,~,~,~,n_st,~,alpha]=parameters(1);
A_prime  = repmat(a',1,na);
A_prime2 = repmat(a',1,na)';
             % /* agents consume their income */
copt = zeros(na,n_st);          % /* optimal consumption */
%aopt = zeros(na,n_st);          % /* optimal next-period assets */

if type==1 || type ==3 %1 and 3 eudated
      endow = endow_all(:,1);
      v	 = zeros(na,n_st); %POS V=NAX2
      v	 = v/(1-beta); 
	  l1 = (endow(1)/(alpha*price))^(1/(alpha-1));
      x1 = max(l1^alpha,0);
      c1 = max(price*x1+(endow(1)*(1-l1))+(1+r).*A_prime2-A_prime,0);
	  
      l2 = (endow(2)/(alpha*price))^(1/(alpha-1));
      x2 = max(l2^alpha,0);
      c2 = max(price*x2+(endow(2)*(1-l2))+(1+r).*A_prime2-A_prime,0);
	  
      l3 = (endow(3)/(alpha*price))^(1/(alpha-1));
      x3 = max(l3^alpha,0);
      c3 = max(price*x3+(endow(3)*(1-l3))+(1+r).*A_prime2-A_prime,0);
	  
	  l4 = (endow(4)/(alpha*price))^(1/(alpha-1));
      x4 = max(l4^alpha,0);
      c4 = max(price*x4+(endow(4)*(1-l4))+(1+r).*A_prime2-A_prime,0);
	  %c1(find(c1<=0)) = NaN;
      %c2(find(c2<=0)) = NaN;
      %c3(find(c3<=0)) = NaN;
      %c4(find(c4<=0)) = NaN;
      util1 = u(c1,[],2);
	  util2 = u(c2,[],2);
	  util3 = u(c3,[],2);
	  util4 = u(c4,[],2);
      %util1(find(isnan(util1))) = -inf;
      %util2(find(isnan(util2))) = -inf;
      %util3(find(isnan(util3))) = -inf;
      %util4(find(isnan(util4))) = -inf;
else
      endow = endow_all(:,2);
      v  = zeros(na,n_st);
      v	 = v/(1-beta); 
	  l1 = (endow(1)/(alpha*price))^(1/(alpha-1));
      x1 = max(l1^alpha,0);
      c1 = max(price*x1+(endow(1)*(1-l1))+(1+r).*A_prime2-A_prime,0);
	  
      l2 = (endow(2)/(alpha*price))^(1/(alpha-1));
      x2 = max(l2^alpha,0);
      c2 = max(price*x2+(endow(2)*(1-l2))+(1+r).*A_prime2-A_prime,0);
	  
      l3 = (endow(3)/(alpha*price))^(1/(alpha-1));
      x3 = max(l3^alpha,0);
      c3 = max(price*x3+(endow(3)*(1-l3))+(1+r).*A_prime2-A_prime,0);
	  
	  l4 = (endow(4)/(alpha*price))^(1/(alpha-1));
      x4 = max(l4^alpha,0);
      c4 = max(price*x4+(endow(4)*(1-l4))+(1+r).*A_prime2-A_prime,0);
      util1 = u(c1,0);
	  util2 = u(c2,0);
	  util3 = u(c3,0);
	  util4 = u(c4,0);

    
end
	  test2 = 1;
	  while (test2 > 10e-10)
         
	  %{
	       r1(:,:)=util1(:,:)+beta*surv(1,1).*repmat((prob(1,1)*v(:,1)+ prob(1,2)*v(:,2)+ prob(1,3)*v(:,3)+ prob(1,4)*v(:,4)),1,nasset);
           r2(:,:)=util2(:,:)+beta*surv(1,2).*repmat((prob(2,1)*v(:,1)+ prob(2,2)*v(:,2)+ prob(2,3)*v(:,3)+ prob(2,4)*v(:,4)),1,nasset);
           r3(:,:)=util3(:,:)+beta*surv(1,3).*repmat((prob(3,1)*v(:,1)+ prob(3,2)*v(:,2)+ prob(3,3)*v(:,3)+ prob(3,4)*v(:,4)),1,nasset);
           r4(:,:)=util4(:,:)+beta*surv(1,4).*repmat((prob(4,1)*v(:,1)+ prob(4,2)*v(:,2)+ prob(4,3)*v(:,3)+ prob(4,4)*v(:,4)),1,nasset);
       %}
   
           r_aux1(:,:)=util1(:,:)+beta.*repmat((pp(1,1)*v(:,1)+ pp(1,2)*v(:,2)+ pp(1,3)*v(:,3)+ pp(1,4)*v(:,4)),1,na);
           r_aux2(:,:)=util2(:,:)+beta.*repmat((pp(2,1)*v(:,1)+ pp(2,2)*v(:,2)+ pp(2,3)*v(:,3)+ pp(2,4)*v(:,4)),1,na);
           r_aux3(:,:)=util3(:,:)+beta.*repmat((pp(3,1)*v(:,1)+ pp(3,2)*v(:,2)+ pp(3,3)*v(:,3)+ pp(3,4)*v(:,4)),1,na);
           r_aux4(:,:)=util4(:,:)+beta.*repmat((pp(4,1)*v(:,1)+ pp(4,2)*v(:,2)+ pp(4,3)*v(:,3)+ pp(4,4)*v(:,4)),1,na);
       

           [tv1,tdecis1]=max(r_aux1);
           [tv2,tdecis2]=max(r_aux2);
           [tv3,tdecis3]=max(r_aux3);
           [tv4,tdecis4]=max(r_aux4);
           tdecis=[tdecis1' tdecis2' tdecis3' tdecis4'];
           tv=[tv1' tv2' tv3' tv4'];

           %test1=max(any(tdecis-decis));
           %test2=max(max(abs(tv - v))');
           test2=max(max(abs((tv-v)./tv)));
           v=tv;
           decis=tdecis;
   end
   aopt=a(decis);
   

		  
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