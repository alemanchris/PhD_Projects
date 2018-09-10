# Only copying sources
 # Aiyagari Model adapted from George Hall code, http://people.brandeis.edu/~ghall/econ303/
     # tested in Julia 0.6
     # this code is part of chapter 6, \ Heterogeneous Agents Models\  from the book: \ Introduction to Quantitative Macroeconomics with Julia\
     # Academic Press - Elsevier
     # for comments, email at: petre(dot)caraiani(at)gmail(dot)com

     #  set parameter values
     sigma  = 1.50;            # risk aversion
     beta   = 0.98;            # subjective discount factor
     prob   = [ .8 .2; .5 .5]; # prob(i,j) = probability (s(t+1)=sj | s(t) = si)
     delta  = 0.97;            # 1 - depreciation
     A      = 1.00;            # production technology
     alpha  = 0.25;            # capital's share of income
     theta  = 0.05;            # non-rental income if unemployed is theta*wage
     Kstart = 10.0;            # initial value for aggregate capital stock
     g      = 0.20;            # relaxation parameter

     #   form capital grid
     maxkap = 20;                      # maximum value of capital grid
     inckap = 0.025;                   # size of capital grid increments
     nkap   = trunc(Int,maxkap/inckap+1);  # number of grid points: make it integer - Julia indexes must be integer

     #global variables
     decis  = zeros(nkap,2);
     lambda = zeros(nkap,2);
     probk  = zeros(nkap,1);

     #   calculate aggregate labor supply
     D = [0.0 0.0; 0.0 0.0];
     ed,ev = eig(prob);

     #   make a matrix from the eigenvalues
     edm   = diagm(ed)
     (emax,inmax) = findmax(edm)

     D[inmax,inmax] = emax;
     pinf = ev*D*inv(ev);
     pempl = pinf[inmax,inmax];
     N = 1.0*pempl + theta*(1-pempl);


    liter   = 1;
     maxiter = 50;
     toler   = 0.001;
     metric  = 10.0;
     K = Kstart;
     Kold= K;
     wage=1.0;
     rent=1.0;
     println( ITERATING ON K  );

     println( Iter    metric                       meanK                        Kold  );


     #   loop to find fixed point for agregate capital stock
     while  (metric[1] > toler) & (liter <= maxiter);


     #  calculate rental rate of capital and wage
        #
        wage = (1-alpha) * A * K^(alpha)   * N^(-alpha);
        rent = (alpha)   * A * K^(alpha-1) * N^(1-alpha);
        #
        #  tabulate the utility function such that for zero or negative
        #  consumption utility remains a large negative number so that
        #  such values will never be chosen as utility maximizing
        #
        util1=-10000*ones(nkap,nkap);  # utility when employed
        util2=-10000*ones(nkap,nkap);  # utility when unemployed
        for i=1:nkap;
              kap=(i-1)*inckap;
              for j=1:nkap;
                    kapp = (j-1)*inckap;
                    cons1 = wage + (rent + delta)*kap - kapp;
                    if cons1 > .0;
                       util1[j,i]=(cons1)^(1-sigma)/(1-sigma);
                    end;
                    cons2 = theta*wage + (rent + delta)*kap - kapp;
                    if cons2 > .0;
                       util2[j,i]=(cons2)^(1-sigma)/(1-sigma);
                    end;
              end;
        end;
        #
        #  initialize some variables
        #
        v       = zeros(nkap,2);
        tdecis1 = zeros(nkap,2);
        tdecis2 = zeros(nkap,2);

        test    = 10;
        rs,cs = size(util1)
        r1=zeros(cs,cs);
        r2=zeros(cs,cs);

        #
        #  iterate on Bellman's equation and get the decision
        #  rules and the value function at the optimum


        while test != 0;
            for i=1:cs;
                r1[:,i]=util1[:,i]+beta*(prob[1,1]*v[:,1]+ prob[1,2]*v[:,2]);
                r2[:,i]=util2[:,i]+beta*(prob[2,1]*v[:,1]+ prob[2,2]*v[:,2]);
            end;
            (tv1,inds1)=findmax(r1,1);
            tdecis1    =map(x->ind2sub(r1, x)[1], inds1)   #to find the relative position of the max in a column
            (tv2,inds2)=findmax(r2,1);
            tdecis2    =map(x->ind2sub(r2, x)[1], inds2)

            tdecis=[tdecis1' tdecis2'];
            tv=[tv1' tv2'];

            test=maximum((tdecis-decis));
            copy!(v, tv);
            copy!(decis, tdecis);
       end;

        decis=(decis-1)*inckap;

     #   form transition matrix
     #   trans is the transition matrix from state at t (row)
     #   to the state at t+1 (column)
     #   The eigenvector associated with the unit eigenvalue
     #   of trans' is  the stationary distribution.

        g2=spzeros(cs,cs);
        g1=spzeros(cs,cs);
        for i=1:cs
            g1[i,tdecis1[i]]=1;
            g2[i,tdecis2[i]]=1;
        end
        trans=[ prob[1,1]*g1 prob[1,2]*g1; prob[2,1]*g2 prob[2,2]*g2];
        trans=trans';
        probst = (1/(2*nkap))*ones(2*nkap,1);
        test=1;
        while test > 10.0^(-8);
           probst1 = trans*probst;
           test = maximum(abs.(probst1-probst));
           copy!(probst, probst1);
        end;

     #   vectorize the decision rule to be conformable with probst
     #   calculate new aggregate capital stock  meanK

     kk=vec(decis);
     meanK=probst'*kk;

     #  calculate measure over (k,s) pairs
     #  lambda has same dimensions as decis
     lambda=reshape(probst, cs,2)


     #   calculate stationary distribution of k
     d1,v1=eig(prob');
     d1m  = diagm(d1)
     dmax,imax=findmax(diag(d1m))

     probst1=v1[:,imax];
     ss=sum(probst1);
     probst1=probst1./ss;
     probk=sum(lambda',1)'
     #   form metric and update K
     Kold= K;
     Knew= g*meanK[1] + (1-g)*Kold;
     metric = abs.((Kold-meanK)./Kold);
     K = Knew;
     println(liter,           ,metric[1],             ,meanK[1],              ,Kold);
     liter = liter+1;
     end;

	  #print  results
     println(  PARAMETER VALUES  );
     println(  sigma      beta      delta       A       alpha      theta  );
     println(sigma,            ,beta,             ,delta,         , A,        ,alpha,       ,theta);
     println(  EQUILIBRIUM RESULTS   );
     println(  K               N              wage             rent  );
     println(round(Kold,4),              , round(N,4),               ,round(wage,4),            , round(rent,4) );

	 # FUNCTIONS

	#  function markov(T,n,s0,V);
    #  chain generates a simulation from a Markov chain of dimension
    #  the size of T
    #  T is transition matrix
    #  n is number of periods to simulate
    #  s0 is initial state
    #  V is the quantity corresponding to each state
    #  state is a matrix recording the number of the realized state at time t

    function markov(T,n,s0,V);

		r,c   = size(T);
		v1,v2 = size(V);

		#rand('uniform');
		#using Distributions
		X=rand(Uniform(0,1), n-1,1)

		state=zeros(2,99);
		chain=[];
		s=zeros(r,1);
		s[s0]=1
		cum=T*triu(ones(size(T)));
		ppi  =[];

		state[:,1]=s

		for k=1:length(X);
		#k=1
		  state[:,k]=s;
		  ppi =[0 s'*cum];
		  ss1= convert(Array{Float64},((X[k].<=ppi[2:r+1])))
		  ss2= convert(Array{Float64},(X[k].>ppi[1:r]))
		  s=(ss1.*ss2)
		  s=reshape(s,2,1)
		end

		chain=V*state
		return state,chain
    end

	# SIMULATION
	 #    simulate life histories of the agent
     #rand('uniform');
     using Distributions

     println(  SIMULATING LIFE HISTORY  );
     k = Kold;               # initial level of capital
     n = 100;                # number of periods to simulate
     s0 = 1;                 # initial state
     hist = zeros(n-1,2);
     cons = zeros(n-1,1);
     invest = zeros(n-1,1);

     grid  = collect(0:inckap:maxkap);
     r,c   = size(prob);
     T     = prob;
     V     = collect(1:r)';

     state, chain= markov(prob,n,s0,V)

     chain=convert(Array{Int64,2},chain);
     state=convert(Array{Int64,2},state);

	 for i = 1:n-1;
         hist[i,:] = [ k chain[i] ];
         I1 = trunc(Int,k/inckap) ;
         I2 = trunc(Int,k/inckap) + 1;
         if I1 == 0;
            I1=1;
            println(  N.B.  I1 = 0  );
         end;
         if I2 > nkap;
            I2 = nkap;
            println(  N.B.  I2 > nkap  );
         end;
         weight = (grid[I2,1] - k)/inckap;
         kprime = weight*(decis[I1,chain[i]]) +  (1-weight)*(decis[I2,chain[i]]);
         if chain[i] == 1;
            cons[i] = wage + (rent + delta)*k - kprime;
         elseif chain[i] == 2;
            cons[i] = wage*theta + (rent + delta)*k - kprime;
         else;
           println(  something is wrong with chain  );
           chain
         end;
         k = kprime;
         invest[i] = kprime;
     end;

	 #plots
     #Income distribution
     income =  [ (rent*grid + wage)  (rent*grid + wage*theta) ]  ;
     #sort income
     pinc = sortrows(hcat(vec(income), 1:length(vec(income))), by = x -> x[1])[:,1]
     index= sortrows(hcat(vec(income), 1:length(vec(income))), by = x -> x[1])[:,2]
     index=trunc.(Int,index)
     plambda = vec(lambda)
	#=
	"using Plots\n",
    "plotly() \n",
    "plot(pinc,plambda[index], linewidth=1,title=\"Income Distribution\", ylabel=\"% of agents\", label=\"income level\")"
	=#
	using Gadfly


    #Capital distribution


    #plot(grid,probk, linewidth=1,title=\"Capital Distribution\", ylabel=\"% of agents\"label=\"capital\")\n"tion with 1 method)
