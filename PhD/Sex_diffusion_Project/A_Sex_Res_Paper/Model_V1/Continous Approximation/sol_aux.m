function H = sol_aux(ax)
H=zeros(2,size(ax,2));
for k=1:size(ax,2)
    r          = ax(1,k);
    price      = ax(2,k);
    rate_price = [r,price];
    excess     = clearing(rate_price);
    
    H(1,k)=excess(1);
    H(2,k)=excess(2);
end


end