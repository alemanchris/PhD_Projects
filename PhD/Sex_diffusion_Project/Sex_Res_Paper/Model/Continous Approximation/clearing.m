function [excess_total,dist]=clearing(rate_price)
% 1 is male
[ex_male_edu,d_m_edu]   = excess(rate_price,1,1); % each d_m has 4 colums
[ex_female_edu,d_f_edu] = excess(rate_price,1,2);
[ex_male_nedu,d_m_nedu]   = excess(rate_price,2,1);
[ex_female_nedu,d_f_nedu] = excess(rate_price,2,2);

% excess demand of assets
ex_as  = ex_male_edu(1)+ex_female_edu(1)+ex_female_nedu(1)+ex_male_nedu(1);
% excess demand of sex
ex_sex = ex_female_nedu(2)-ex_male_nedu(2)+ex_female_edu(2)-ex_male_edu(2);
% Excess vector
excess_total = [ex_as,ex_sex];
if nargout>1
    dist = [d_m_edu,d_f_edu,d_m_nedu,d_f_nedu]; %16 column
    dist = dist/sum(sum(dist));
else 
end

end