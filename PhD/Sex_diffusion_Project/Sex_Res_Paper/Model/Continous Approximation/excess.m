function [ex_sa,distrib]=excess(rate_price,edu,gender)
[~,~,~,a,na,~,pp,~,~,pp1,amin1,amax1,~,nk,~,ag,n_st,ngk,~,~,~,~,~,gamma_all,phy]=parameters(1);
r  = rate_price(1);
price = rate_price(2); 
if edu ==1
    gamma = gamma_all(:,1);
else
    gamma = gamma_all(:,2);
end
if gender==1
[aopt,~,xopt,~]=partial_gs(r,price,edu,gender);
else
[aopt,~,xopt,~]=partial_gs(r,price,edu,gender);
end
     q1    = 0;   % Iteration counter for invariant distribution 
	  kritg = 1;   % Remember we will be updating kritg as kritg = sum(abs(gk0-gk));
	  % Populate the distributions
	  % Initialization of the distribution functions
	  %if q<=10 % q is the interest rate iteration counter
		gk = ones(nk,n_st)/nk; %NK x 2
        % This makes a 300X2.*1X2
		gk = gk.*pp1'; %POS pp1:2x1 is the invariant of the transition matrix
        % does it really matter to multiply by pp1? does it speed things
        % up?
	  %end
	  %if q==10
	%	ngk = 5*ngk;
	 % end
	  %disp ()
	   % == es <
	 % > es <a+1 or <=
	 % < es >= or >a-1
      while kritg>10e-10
	  %while  q1<=ngk
			q1 = q1+1;
			gk0 = gk; %Remember gk is a 300X2 that was multiplied by the ergodic
			gk  = zeros(nk,n_st); % reiniate gk
			% l=i
			% i=ii
			for i=1:n_st
				for ii=1:nk %nk is the number of elements on AG
					k0 = ag(ii); % Running element of ag
					if k0<=amin1 % if im on the the left limmit which happens always with equality
						k1 = aopt(1,i);
					elseif k0>=amax1 % if im on the right limit which happens always with equality
						k1 = aopt(na,i);
					else 
						k1 = interp1(a,aopt(:,i),k0);
						%k1 = lininterp(a,aopt(:,i),k0,'spline');
					end
					
					if k1<=amin1 
                        % as a reference there is the pp= [0.925 0.075;0.5 0.5];
                        % pp(i,1) if i=1 is the prob of being employed and keep ur
                        % job 
                        % pp(i,2) if i=1 if the prob of employ and loosing.
						% pp(2,2) if i=1 if the prob of unemploy and staying unemploy.
						gk(1,1) = gk(1,1) + gk0(ii,i)*pp(i,1)+(0.5*phy(edu,1)*gk0(1,1));
						gk(1,2) = gk(1,2) + gk0(ii,i)*pp(i,2)+(0.5*phy(edu,2)*gk0(1,2));
                        gk(1,3) = gk(1,3) + gk0(ii,i)*pp(i,3)+(0.5*phy(edu,3)*gk0(1,3));
						gk(1,4) = gk(1,4) + gk0(ii,i)*pp(i,4)+(0.5*phy(edu,4)*gk0(1,4));
                        % Esto llena las distribuciones fila por fila. (pero solo la primera)
					elseif k1>=amax1
						gk(nk,1) = gk(nk,1) + gk0(ii,i)*pp(i,1)+(0.5*phy(edu,1)*gk0(nk,1));
						gk(nk,2) = gk(nk,2) + gk0(ii,i)*pp(i,2)+(0.5*phy(edu,2)*gk0(nk,2));
                        gk(nk,3) = gk(nk,3) + gk0(ii,i)*pp(i,3)+(0.5*phy(edu,3)*gk0(nk,3));
						gk(nk,4) = gk(nk,4) + gk0(ii,i)*pp(i,4)+(0.5*phy(edu,4)*gk0(nk,4));
					elseif (k1>amin1) && (k1<amax1)
						%j = sum(ag'.<=k1)+1; %in GAUSS
                        % ag is 1X81, so j tells me the number of points in
                        % the grid that are equal or below k1+1
                        % it is +1 because it is a position counter to
                        % start filling
                        j = sum(ag'<=k1)+1; % 
                        % keep in minf that k1> ag(j-1)always
						n0 = (k1-ag(j-1))/(ag(j)-ag(j-1)); % Esta es la clave, de por que le sale bonita la distribucion, por que cada bin lo divide en dos bins.
                        % The double gk is a way to sum cuz in the first
                        % loop it stores in gk(j,1) the value of the first
                        % loop over i(employed)
						% Es como decir : 
                        % n0*gk0(ii,1)*pp(1,1)+n0*gk0(ii,2)*pp(2,1)
						gk(j,1) = gk(j,1)+n0*gamma(i)*gk0(ii,i)*pp(i,1)+(0.5*phy(edu,1)*gk0(j,1)); %starts filling k1+1 (I multiply by 0.5 cuz it double sums so at the end I only sum once)
						gk(j,2) = gk(j,2)+n0*gamma(i)*gk0(ii,i)*pp(i,2)+(0.5*phy(edu,2)*gk0(j,2));
                        gk(j,3) = gk(j,3)+n0*gamma(i)*gk0(ii,i)*pp(i,3)+(0.5*phy(edu,3)*gk0(j,3)); %starts filling k1+1
						gk(j,4) = gk(j,4)+n0*gamma(i)*gk0(ii,i)*pp(i,4)+(0.5*phy(edu,4)*gk0(j,4));
						gk(j-1,1) = gk(j-1,1)+(1-n0)*gamma(i)*gk0(ii,i)*pp(i,1)+(0.5*phy(edu,1)*gk0(j,1));
						gk(j-1,2) = gk(j-1,2)+(1-n0)*gamma(i)*gk0(ii,i)*pp(i,2)+(0.5*phy(edu,2)*gk0(j,2));
                        gk(j-1,3) = gk(j-1,3)+(1-n0)*gamma(i)*gk0(ii,i)*pp(i,3)+(0.5*phy(edu,3)*gk0(j,3));
						gk(j-1,4) = gk(j-1,4)+(1-n0)*gamma(i)*gk0(ii,i)*pp(i,4)+(0.5*phy(edu,4)*gk0(j,4));
						% add new born

					end
				end
			end
			
			gk 	  = gk/sum(sum(gk));
            %kritg = sum(abs(gk0-gk));
			kritg = max(max(abs(gk0-gk)));
	  
	 end
	 a0 = (gk(:,1)+gk(:,2)+gk(:,3)+gk(:,4))'*ag';   % This is the average asset level
     s0 = (interp1(a',sum(xopt,2),ag')')*sum(gk,2); % Excess sex consumption 
     %Vq = interp1(X,V,Xq) interpolates to find Vq, the values of the
     %underlying function V=F(X) at the query points Xq. 
     ex_s = s0;
     ex_a = a0;
     ex_sa=[ex_a,ex_s];
     if nargout>1
         distrib=gk;
     else
     end
end

%{
      q1    = 0;   % Iteration counter for invariant distribution 
	  kritg = 1;   % Rember we will updating kritg as kritg = sum(abs(gk0-gk));
	  % Populate the distributions
	  % Initialization of the distribution functions
	  %if q<=10 % q is the interest rate iteration counter
		gk = ones(nk,n_st)/nk; %NK x 2
        % This makes a 300X2.*1X2
		gk = gk.*pp1'; %POS pp1:2x1 is the invariant of the transition matrix
        % does it really matter to multiply by pp1? does it speed things
        % up?
	  %end
	  %if q==10
		%ngk = 5*ngk;
	  %end
	  %disp ()
	   % == es <
	 % > es <a+1 or <=
	 % < es >= or >a-1
      %while kritg>10e-10
	  while  q1<=ngk
			q1 = q1+1;
			gk0 = gk; %Remember gk is a 300X2 that was multiplied by the ergodic
			gk  = zeros(nk,n_st); % reiniate gk
			% l=i
			% i=ii
			for i=1:n_st
				for ii=1:nk %nk is the number of elements on AG
					k0 = ag(ii); % Running element of ag
					if k0<=amin1 % if im on the the left limmit which happens always with equality
						k1 = aopt(1,i);
					elseif k0>=amax1 % if im on the right limit which happens always with equality
						k1 = aopt(na,i);
					else 
						%k1 = interp1(,1)
						k1 = lininterp(a,aopt(:,i),k0);
					end
					
					if k1<=amin1 
                        % as a reference there is the pp= [0.925 0.075;0.5 0.5];
                        % pp(i,1) if i=1 is the prob of being employed and keep ur
                        % job 
                        % pp(i,2) if i=1 if the prob of employ and loosing.
						gk(1,1) = gk(1,1) + gk0(ii,i)*pp(i,1);
						gk(1,2) = gk(1,2) + gk0(ii,i)*pp(i,2);
                        gk(1,3) = gk(1,3) + gk0(ii,i)*pp(i,3);
						gk(1,4) = gk(1,4) + gk0(ii,i)*pp(i,4);
                        % Esto llena las distribuciones fila por fila. (pero solo la primera)
					elseif k1>=amax1
						gk(nk,1) = gk(nk,1) + gk0(ii,i)*pp(i,1);
						gk(nk,2) = gk(nk,2) + gk0(ii,i)*pp(i,2);
                        gk(nk,3) = gk(nk,3) + gk0(ii,i)*pp(i,3);
						gk(nk,4) = gk(nk,4) + gk0(ii,i)*pp(i,4);
					elseif (k1>amin1) && (k1<amax1)
						%j = sum(ag'.<=k1)+1; %in GAUSS
                        % ag is 1X81, so j tells me the number of points in
                        % the grid that are equal or below k1+1
                        % it is +1 because it is a position counter to
                        % start filling
                        j = sum(ag'<=k1)+1; % 
                        % keep in minf that k1> ag(j-1)always
						n0 = (k1-ag(j-1))/(ag(j)-ag(j-1));
                        % The double gk is a way to sum cuz in the first
                        % loop it stores in gk(j,1) the value of the first
                        % loop over i(employed)
                        % I still suspect it is double doing per j
						gk(j,1) = gk(j,1)+n0*gk0(ii,i)*pp(i,1); %starts filling k1+1
						gk(j,2) = gk(j,2)+n0*gk0(ii,i)*pp(i,2);
                        gk(j,3) = gk(j,3)+n0*gk0(ii,i)*pp(i,3); %starts filling k1+1
						gk(j,4) = gk(j,4)+n0*gk0(ii,i)*pp(i,4);
						gk(j-1,1) = gk(j-1,1)+(1-n0)*gk0(ii,i)*pp(i,1);
						gk(j-1,2) = gk(j-1,2)+(1-n0)*gk0(ii,i)*pp(i,2);
                        gk(j-1,3) = gk(j-1,3)+(1-n0)*gk0(ii,i)*pp(i,3);
						gk(j-1,4) = gk(j-1,4)+(1-n0)*gk0(ii,i)*pp(i,4);
					end
				end
			end
			gk 	  = gk/sum(sum(gk));
            %kritg = sum(abs(gk0-gk));
			kritg = max(max(abs(gk0-gk)));
	  
	 end
	 a0 = (gk(:,1)+gk(:,2)+gk(:,3)+gk(:,4))'*ag';  % This is the average asset level
     s0 = (interp1(a',sum(xopt,2),ag')')*sum(gk,2); % Excess sex consumption 
     %Vq = interp1(X,V,Xq) interpolates to find Vq, the values of the
     %underlying function V=F(X) at the query points Xq. 
     ex_s = s0;
     ex_a = a0;
     ex_sa=[ex_a,ex_s];
     if nargout>1
         distrib=gk;
     else
     end
end
%}