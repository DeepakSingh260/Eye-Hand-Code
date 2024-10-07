% sampDC.m
%
% usage: Dcell=sampDC(Dcell);

function Dcell=sampDC(Dcell)
global Sexp Etmp DCtmp
if nargin==1, sampdc=Dcell(end,1:2); Dcell(end,:)=[]; else sampdc=Etmp; end
for Enow=Sexp.EYEnow, sampdc{Enow}=sampdc{Enow}(1:Sexp.jE,:); end
Etmp=cell(1,2); Sexp.Ccoords{3}(Sexp.Lc)=GetSecs; %clear leftover global dataset and timestamp the new drift correction

Sexp.DCdat(Sexp.Lc,:)=cell(1,2);
for Enow=Sexp.EYEnow, 
	if or(nargin==0,Sexp.Lc==1), lastcorrect=[0 0]; else lastcorrect=Sexp.Ccoords{Enow}(Sexp.Lc-1,1:2); end %extract last correction
	sampdc{Enow}=sampdc{Enow}(1:find(~isnan(sampdc{Enow}(:,1)),1,'last'),:); %extract zeroed&flipped central fixation data
	Sexp.DCdat{Sexp.Lc,Enow}=[sampdc{Enow}(:,1)+lastcorrect(1) sampdc{Enow}(:,2)-lastcorrect(2) sampdc{Enow}(:,3:4)]; %re-introduce last correction
	Sexp.Ccoords{Enow}(Sexp.Lc,:)=nanmedian(Sexp.DCdat{Sexp.Lc,Enow}(find(Sexp.DCdat{Sexp.Lc,Enow}(:,4)<=Sexp.DCdat{Sexp.Lc,Enow}(Sexp.jE,4),10,'last'),1:2),1); end %new Ccoords is median of (zero-centered) fixation positions (pos-right, pos-up)
Sexp.Lc=Sexp.Lc+1; DCtmp=[];
