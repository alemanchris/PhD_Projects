function [util]=u(x,y,gender)
[sigma]=parameters(1);
if gender==2 % Producers
util = (x.^(1-sigma))./(1-sigma);
else
util = (x.^(1-sigma))./(1-sigma)+(y.^(1-sigma))./(1-sigma);
end
end