%windowmean.m
%
% [meanvec]=windowmean(dat,windowsize,FigNum,varargin);
%
% dat:        vector of length N
% windowsize: size of sliding window for mean calc
% Fignum:     Leave empty for no plot. 0: plot new fig. NonZero: plot figure(plotNY01)
% COL:        Color: either a 2-row, 3-column matrix of (row1) markerface
%                    (row2) markeredge colors, or two characters for face
%                    and edge colors. Can also specify COLedge and COLface
%                    [note: for the dot symbol, the color is specified as
%                    an edge, not a face color]
% causal:     Causal filter, or centered window (with nan-padding)
% skip:       Return data indices at window centers for non-overlapping windows.
%             Centers the set of returned datapoints on the dataset, while
%             leaving equal (to within one) nan-padded estinates at both
%             ends. Produces about length(dat)/windowsize datapoints.
%             skip=[][no skip - default]; 
%             skip=0[center the skipped trace (maximum symmetrical coverage)]; 
%             skip=1[first mean at t=1];
%             skip=2[last mean at t=length(dat)];
%
% teh wrote it. Oct 11, 2009

% You can use filter to find a running average without using a for loop. This example finds the running average of a 16-element vector, using a window size of 5.
% 
% data = [1:0.2:4]';
% windowSize = 5;
% filter(ones(1,windowSize)/windowSize,1,data)


function [meanvec inds]=windowmean(dat,windowsize,varargin)
if nargin==1, windowsize=round(length(dat)/40); end

va=varargin; FigNum=0; COLedge=[]; COLface=[]; Causal=[]; Skip=[];
if ~isempty(va), if isnumeric(va{1}), FigNum=va{1}; va=va(2:end); else
	for n=1:length(va), if ischar(va{n}), if and(or(strcmp(va{n},'FigNum'),strcmp(va{n},'fignum')),length(va)>n), if or(isempty(va{n+1}),and(isnumeric(va{n+1}),all(size(va{n+1})==1))), FigNum=va{n+1}; else error('Incorrect assignment to ''FigNum'''); end, end, end, end, end
for n=1:length(va), if ischar(va{n}), if and(or(strcmp(va{n},'COL'),strcmp(va{n},'col')),length(va)>n), if and(ischar(va{n+1}),length(va{n+1})<=2), COLface=va{n+1}(1); if length(va{n+1})==2, COLedge=va{n+1}(2); else COLedge='k'; end, elseif and(isnumeric(va{n+1}),size(va{n+1},2)==3), COLface=va{n+1}(1,:); if size(va{n+1},1)==2, COLedge=va{n+1}(2,:); else COLedge='k'; end, else error('Incorrect assignment to property ''COL'''); end, end, end, end
for n=1:length(va), if ischar(va{n}), if and(or(strcmp(va{n},'COLedge'),strcmp(va{n},'coledge')),length(va)>n), if and(ischar(va{n+1}),length(va{n+1})==1), COLedge=va{n+1}; elseif and(isnumeric(va{n+1}),all(size(va{n+1})==1)), COLedge=va{n+1}*[1 1 1]; elseif and(isnumeric(va{n+1}),and(size(va{n+1},1)==1,size(va{n+1},2)==3)), COLedge=va{n+1}; else error('Incorrect assignment to property ''COLedge'''); end, end, end, end
for n=1:length(va), if ischar(va{n}), if and(or(strcmp(va{n},'COLface'),strcmp(va{n},'colface')),length(va)>n), if and(ischar(va{n+1}),length(va{n+1})==1), COLface=va{n+1}; elseif and(isnumeric(va{n+1}),all(size(va{n+1})==1)), COLface=va{n+1}*[1 1 1]; elseif and(isnumeric(va{n+1}),and(size(va{n+1},1)==1,size(va{n+1},2)==3)), COLface=va{n+1}; else error('Incorrect assignment to property ''COLface'''); end, end, end, end
for n=1:length(va), if ischar(va{n}), if and(or(strcmp(va{n},'Causal'),strcmp(va{n},'causal')),length(va)>n), if and(isnumeric(va{n+1}),length(va{n+1})==1), Causal=va{n+1}; else error('Incorrect assignment to property ''Causal'''); end, end, end, end
for n=1:length(va), if ischar(va{n}), if and(or(strcmp(va{n},'Skip'),strcmp(va{n},'skip')),length(va)>=n+1), if and(isnumeric(va{n+1}),all(size(va{n+1})==1)), Skip=va{n+1}; if ~or(Skip==0,or(Skip==1,Skip==2)) error('Incorrect assignment to ''FigNum'''); end, end, end, end, end, end
if isempty(COLedge), COLedge='k'; end
if isempty(COLface), COLface=[.5 .5 .5]; end
if isempty(Causal), Causal=0; end

if Causal, nanlength=windowsize-1; else nanlength=floor(windowsize/2); end
keeperinds=1:length(dat);
dat=[nan(nanlength,1); dat(:)];
datmat=nan(length(dat),windowsize);

for c=1:windowsize, datmat(1:length(dat),c)=dat; dat=dat(2:end); end
DEN=sum(~isnan(datmat),2); datmat(isnan(datmat))=0; meanvec=sum(datmat,2)./DEN; meanvec=meanvec(keeperinds);

inds=[1:length(meanvec)]'; 
if Skip>0, inds=Skip:windowsize:length(meanvec); 
elseif Skip==0, inds=1:windowsize:length(meanvec); di=length(meanvec)-inds(end); inds=inds+floor(di/2); end
if ~isempty(FigNum),
	if FigNum==0, figure; else figure(FigNum); hold on; end
	if length(meanvec)>100, plot(inds,meanvec(inds),'.','MarkerEdgeColor',COLedge); 
	else plot(inds,meanvec(inds),'ko--','MarkerFaceColor',COLface,'MarkerEdgeColor',COLedge,'MarkerSize',7); end, end
