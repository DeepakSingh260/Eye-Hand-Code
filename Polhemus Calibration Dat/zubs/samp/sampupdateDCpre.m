% sampupdateDCpre.m
% update Etmp such that it maintains same length, while
%   simultaneously testing whether all data are within Dthresh of
%   datamedian
%
% usage: goodtogo=sampupdateDCpre(Dthresh);

function goodtogo=sampupdateDCpre(Dthresh)
global Sexp Etmp DCtmp

Ltrace=111;
if Sexp.jE>Ltrace,
    Sexp.jE=Sexp.jE-1;
    if isempty(DCtmp),
        DCtmp=nan(Sexp.jE,2,2);
        for Enow=Sexp.EYEnow, DCtmp(:,:,Enow)=Etmp{Enow}(1:Sexp.jE,1:2); end
        DCtmp=nanmean(DCtmp,3); end
    
    DCtmpnew=nan(2);
    for Enow=Sexp.EYEnow,
        Etmp{Enow}(1:Sexp.jE,:)=Etmp{Enow}(2:Sexp.jE+1,:);
        DCtmpnew(Enow,:)=Etmp{Enow}(Sexp.jE,1:2); end
    DCtmp=[DCtmp(2:end,:); nanmean(DCtmpnew,1)];
    goodtogo=sum(dist(DCtmp,median(DCtmp,1),2)>Dthresh)<.15*Ltrace; %allow for bet 10-15/100 outliers
else goodtogo=0; end