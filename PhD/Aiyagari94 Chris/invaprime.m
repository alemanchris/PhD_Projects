function [resid]=invaprime(a,y,r,a0,cp0,aprime,sy,w0)
      %aprime = interp1(a0(:,sy),Amat(:,sy),a,'pchip');
      cpol   = interp1(a0(:,sy),cp0(:,sy),a,'pchip');
      resid = (1+r)*a+y(sy)*w0-cpol-aprime;
end