% pollEL.m
%
% usage: [dataavalable]=pollEL; 
%            dataavalable=1 if there is a data update
%            dataavalable=2 if the update is non-nan

function dataavalable=pollEL
global Sexp Etmp ELxy %Etmp is the data, now as a 2-element cell [x,y,diam,t,t], ELxy=[x,y,t,v]
t=GetSecs;

% get data in the form of an event structure
%	 1: time of sample (when camera imaged eye, in milliseconds since tracker was activated)
%   12: left pupil size (arbitrary units
%   13: right pupil size (arbitrary units
%   14: left gaze position x (in pixel coordinates set by screen_pixel_coords command)
%   15: right gaze position x (pix)
%   16: left gaze position y (pix)
%   17: right gaze position y (pix)
%	18: angular resolution x (at gaze position in screen pixels per visual degree)
%	19: angular resolution y (at gaze position in screen pixels per visual degree)

loops=1; tmpcell=cell(2,1); drained=1;
while drained==1, 
    [tmpcell{loops},~,drained]=Eyelink('GetQueuedData');
	if ~isempty(tmpcell{loops}),
		tmpcell{loops}=tmpcell{loops}(:,tmpcell{loops}(1,:)~=-32768)';
		loops=loops+1; end, end
%neeed to cell2mat the cell vect, and then sort using:
newmat=cell2mat(tmpcell);
if ~isempty(newmat),
    newmat(newmat==-32768)=nan;
    newmat=sortrows(newmat,1);
    newmat(all(isnan(newmat(:,14:15)),2),:)=[]; end

sz=size(newmat,1); 
if sz>0, dataavalable=2; ELxy=nan(1,2); Sexp.jE=Sexp.jE+[1:sz]; 
    OST=t+(newmat(:,1)-max(newmat(:,1)))/1000; %OST uses time relative in sec to correct when multiple data per poll
	for Enow=Sexp.EYEnow,
		Etmp{Enow}(Sexp.jE,1:Sexp.nsampEL-1)=[newmat(:,Sexp.LRElist(:,Enow)) OST]; %X/Y/PA/T/v [note vel is computed only at trial end]
		ELxy=nanmean([ELxy; Etmp{Enow}(Sexp.jE(end),[1 2])],1); end
	ELxy=[ELxy t]; Sexp.jE=Sexp.jE(end);
else dataavalable=0; end

