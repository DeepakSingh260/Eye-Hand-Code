%RFL.m remove from list
% Remove a specified set of ROW elements from the input matrix
%   Mout: shortened matrix;
%   Mrem: removed matrix elements
%
% Usage: [Mout Mremoved]=RFL(Min,listtoremove)

function [Mout Mremoved]=RFL(Min,list,keeporient)
if nargin<3, keeporient=0; end
if and(~keeporient,size(Min,1)==1), Min=Min(:); end

%new, faster version
try iis=ones(size(Min,1),1); iis(list)=0;
	Mout=Min(logical(iis),:,:);
	if nargout==2, Mremoved=Min(~logical(iis),:,:); end
catch %use old (slower) RFL if thee is an error
	Mout=Min; list=sort(list);
	Mremoved=Mout(list,:);
	for n=length(list):-1:1
		Mout=[Mout(1:list(n)-1,:,:); Mout(list(n)+1:end,:,:)]; end, end, end