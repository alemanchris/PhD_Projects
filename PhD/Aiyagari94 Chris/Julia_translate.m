function invariant = compute_invariant_pwlinear(N,agrid_finer,Y,agrid,endog_grid,K,L_invariant) 
disp('calculating the invariant distribution via piecewise linear interpolation')
%ay.ny = N
%ay.ny = agrid_finer
%ay.agrid
%ay.a_ast_mat % Endogenous grid
%ay.ygrid
  Ln_mat   = zeros(K, N) ;
  Lnm1_mat = L_invariant; 

  iter = 0; 
  dist = 10; 
 
  
  
  %{
  
countries = ("Japan", "Korea", "China")
cities = ("Tokyo", "Seoul", "Beijing")
for (i, country) in enumerate(countries)
    city = cities[i]
    println("The capital of $country is $city")
end

  %}
%a_ast_itp = interpolate((agrid,Y),endog_grid,Gridded(Linear()));% # next period's assets, current y 
%max_value = max(max(endog_grid));
a_ast_itp = @(a,b) interp2(agrid,Y,endog_grid,a,b,'spline'); % whi not linear, cuz NaN will always be ataken as a minimum
 
    while iter<2000 && dist>1e-4 
 
       iter = iter+1 

 
      for i_y = 1:N %# next period 
           for i_a =1:K %# next period 
              % a_v is the value
                    % i_a is the index
              a_v = agrid_finer(i_a);
              vec_temp = zeros(N,1); 

                 for i_y0 = 1:N   % # last period y 
                     y0_v = Y(i_y0);
                     aval  = min(max(a_ast_itp(a_v, y0_v), agrid(1)), agrid(end));% # today's assets (endogenous grid) 
                     %ind_r = min(max(searchsortedfirst(ay.agrid_finer, aval), 2), K) 
                     ind_r = min(max(nearestpoint(aval,agrid_finer,'nearest'),2),K); % nice trick
 
                     Lval = Lnm1_mat(ind_r-1,i_y0) + (Lnm1_mat(ind_r, i_y0)- Lnm1_mat(ind_r-1, i_y0)) / (agrid_finer(ind_r) - agrid_finer(ind_r-1)) * (aval - agrid_finer(ind_r-1)); 
                     vec_temp(i_y0) = Lval; 
                 end 
 
 
                 %Ln_mat(i_a, i_y) = dot(pi(:, i_y),vec_temp); 
                 Ln_mat(i_a, i_y) = dot(pi(:, i_y),vec_temp); 
             end 
         end 

         dist = maxabs(Ln_mat - Lnm1_mat); 
         Lnm1_mat=Ln_mat;
         %copy!(?nm1_mat, ?n_mat) 

         %# if iter%200 == 0.0 
         %#     @printf("Iteration # %d on distribution with distance %.3g\n", iter, dist) 
         %# end 
    end 
 
    invariant = Lnm1_mat;
    %copy!(?_invariant, ?nm1_mat) 
 
 
   % Void 
 end 
