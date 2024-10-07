%nearest.m
% Finds the index of the nearest value in list to the value of comp
%
% USAGE: ind=nearest(list,comp)
% 
% INPUTS:
% list: either a vector or matrix of values to be compared to comp
% comp: value that is compared to every element in list to determine the closest match
%       NOTE: comp may have multiple values; in this case it should be arranged in a column
%
% OUTPUT:
% ind: the index or indices for the nearest value to comp found in list. Note that if
%      there are multiple nearest-matches and more than one comp value to be tested
%      ind will be a cell vector 
%
% teh wrote it. 3.14.15

function indOut=findnearestN(list,comp,N) %#ok<*TRYNC>
if and(size(list,1)==1,size(list,2)>1), list=list(:); end
if size(list,2)==1, list=[list(:) zeros(length(list),1)]; comp=[comp 0]; end
if nargin==2, N=3; end

indOut=[];
dif=dist(list,comp,2);
while length(indOut)<N,
	dif(indOut)=max(dif)+1;
	indOut=[indOut; str8n(find(isnear(dif,min(dif(:)))))]; end
indOut=sort(indOut(1:N));