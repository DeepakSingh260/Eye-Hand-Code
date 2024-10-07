function h=plotmat(MATin,varargin) %#ok<STOUT,*AGROW>

strarg=[];
for n=1:length(varargin),
	if ischar(varargin{n}), strarg=[strarg '''' varargin{n} ''',']; 
	elseif length(varargin{n})>1, strarg=[strarg '[' num2str(varargin{n}) '],'];
	else strarg=[strarg num2str(varargin{n}) ',']; end, end
for N=1:size(MATin,3),
	if size(MATin,2)==2,
		eval(['h=plot(MATin(:,1,N),MATin(:,2,N),' strarg(1:end-1) ');']);
	elseif size(MATin,2)==3,
		eval(['h=plot3(MATin(:,1,N),MATin(:,2,N),MATin(:,3,N),' strarg(1:end-1) ');']); end, end
if nargout==0, clear h; end