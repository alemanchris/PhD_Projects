function [IND, D] = nearestpoint(x, y, m)
%{
IND contains the indices of each of these points. 
Example: 
 NEARESTPOINT([1 4 12],[0 3]) -> [1 2 2] 
%}
narginchk(2,3);
if nargin ==2
	m = 'nearest';
else 
	if ~ischar(m)
		error('Some Error');
	end
end

if ~isa(x,'double')||~isa(y,'double')
	error('should be double');
end
if isempty(x)||isempty(y)
	IND=[];
	D = [];
	return ;
end

% sort the imput vectors
sz = size(x);
[x,xi]=sort(x(:));
[~,xi]=sort(xi);
nx = numel(x);
cx=zeros(nx,1);
qx = isnan(x);
[y,yi]=sort(y(:));
ny = length(y);
cy = ones(ny,1);
xy = [x;y];
[~,xyi]=sort(xy);
cxy=[cx;cy];
cxy = cxy(xyi);
ii=cumsum(cxy);
ii = ii(cxy==0).';

clear cxy xy xyi;
switch lower(m)
	case {'nearest','near','absolute'}
		% the indices of the nearest point
		ii = [ii;ii+1];
		ii(ii==0)=1;
		ii(ii>ny)=ny;
		yy=y(ii);
		dy = abs(repmat(x.',2,1)-yy);
		[~,ai] = min(dy);
		IND = ii(sub2ind(size(ii),ai,1:nx));
	case {'previous','prev','before'}
		ii(ii<1)=NaN;
		IND =ii;
	case{'next','after'}
		ii=ii+1;
		ii(ii>ny)=NaN;
		IND = ii;
	otherwise
		error('unknown method');
end

IND(qx)=NaN;
if nargout ==2
	D=NaN(1,nx);
	q=~isnan(IND);
	if any(q)
		D(q)=abs(x(q)-reshape(y(IND(q)),[],1));
	end
	D= reshape(D(xi),sz);
end
IND= reshape(IND(xi),sz);
q=~isnan(IND);
IND(q)=yi(IND(q));

end