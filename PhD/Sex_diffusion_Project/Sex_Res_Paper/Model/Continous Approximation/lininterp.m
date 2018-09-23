function [int]=lininterp(xd,yd,x)
% xd = a'
% yd = vold
% x = a1 = 14;
 j = sum(xd<=x');
 %j = sum(xd.<=x');
 int = yd(j)+(yd(j+1)-yd(j)).*(x-xd(j))./(xd(j+1)-xd(j));

end