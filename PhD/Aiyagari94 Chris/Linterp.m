function [Lam] = Linterp(k,i,L0,a,agrid2)

    Lam =  L0(k,i)+(L0(k+1,i)-L0(k,i))*(a-agrid2(k))/(agrid2(k+1)-agrid2(k));

end