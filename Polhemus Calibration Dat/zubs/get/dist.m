% dist.m
% usage:
% [d]=dist(v1,v2,dim)
%
% v1: first input vector/matrix - arranged so that size(v1,dim)=2 for 2 dimensional distance
% v2: second input vector/matrix - arranged so that size(v2,dim)=2 for 2 dimensional distance
% dim: The DIM input identifies the dimension along which the the sum over
%      the squared coordinate elements is taken (to compute D2)
% d:  distance output - ALWAYS ARRANGED IN A COLUMN
%
% teh wrote it (Oct 2, 2009)

function d=dist(v1,v2,dim)
d=[];
if nargin==3,
	if size(dim,2)==2, 
		L=[size(v1,1)==1 size(v2,1)==1];
		if and(L(1),~L(2)), v1=ones(size(v2,1),1)*v1;
		elseif and(L(2),~L(1)), v2=ones(size(v1,1),1)*v2; end
		
		dirvec=dim/sqrt(dot(dim,dim));
		d=dot(v2-v1,ones(size(v2,1),1)*dirvec,2); end, end
if isempty(d),
	if nargin==2,
		if all(size(v2)==1), dim=v2; v2=zeros(size(v1));
		else if all(size(v1)==size(v2')), v2=v2'; end
			if size(v1,2)==2, dim=2; else dim=1; end, end
	elseif nargin==1, v2=zeros(size(v1));
		if size(v1,2)==2, dim=2; else dim=1; end,
	elseif and(nargin==3,and(size(v1,1)>1,size(v2,1)==1)), v2=ones(size(v1,1),1)*v2;
	elseif and(nargin==3,and(size(v1,1)==1,size(v2,1)>1)), v1=ones(size(v2,1),1)*v1; end
	
	if and(size(v1,1)>1,size(v2,1)==1), v2=ones(size(v1,1),1)*v2; end
	
	d=sqrt(sum((v1-v2).^2,dim)); d=d(:); end
