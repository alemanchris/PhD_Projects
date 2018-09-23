function [excess_total,dist]=clearing(rate_price)
% 1 is male
[ex_male,d_m]   = excess(rate_price,1,1);
[ex_female,d_f] = excess(rate_price,1,2);
% excess demand of assets
ex_as  = ex_male(1)+ex_female(1);
% excess demand of sex
ex_sex = ex_female(2)-ex_male(2);
% Excess vector
excess_total = [ex_as,ex_sex];
if nargout>1
    dist = [d_m,d_f];
    dist = dist/sum(sum(dist));
else 
end

end