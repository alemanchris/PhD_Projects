function [inv_ec] = equivec1(p)
 [n,c]=size(p);
 % The aux lines, extract the diagonal of a matrix and plug a different one
 aux  = diag(p)-ones(n,1);
 aux2 = p-diag(diag(p),0); % plug diag(p) on the main diagonal of a matrix of zeros matrix(thats why its zero)
 aux3 = diag(aux,0);
 p =  aux2+aux3;
 p = [p(:,1:n-1) ones(n,1)];
 x = [zeros(n-1,1);1];
 inv_ec = (x'*inv(p))';
end