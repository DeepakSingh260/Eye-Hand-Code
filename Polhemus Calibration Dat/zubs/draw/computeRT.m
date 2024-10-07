%computeRT.m
% used to determine the target radius that will produce a given hit rate
% based on practice data [NB: superior to rT=raylfit(dist(DAT,2),DesiredHitRate)]
%
% usage: rT=computeRT(DAT,DesiredHitRate);
%
% Inputs
%            DAT: 2D dataset 
%                 input as DAT-mean(DAT,1) to compute rel to endpt precision
%                 input as DAT-target to compute rel to accuracy+precision
% DesiredHitRate: desired proportion of endpoints falling within rT
%           Rmax: (optional) maximum allowable target radius
%           Fnum: (optional) produces a plot of the result at figure(Fnum)
%
% Example:
% DAT=mvnrnd([0 0],[14 0;0 14],250);
% DesiredHitRate=.62;
% rT=computeRT(DAT,DesiredHitRate,[],1)
% 
% teh wrote it [10.17.12]

function [rT rTrange]=computeRT(DAT,DesiredHitRate,Rmax,Fnum)
%NB: rTrange will be equal to TOL unless Ncrit terminates the search
if nargin<3, Rmax=[]; end
if nargin<4, Fnum=[]; end

TOL=.001; Ncrit=100;
Dlist=sort(dist(DAT,2));
Robs=[0:length(Dlist)-1]'/length(Dlist);

%FIND correct range of Bhat values
Bhat=.005*[.5 1]; SSD=1;
while SSD>0,
	Bhat=2*Bhat;
	Rhat=raylcdf(Dlist,Bhat(2));
	SSD=sum(Rhat-Robs); end

%PLOT the results from top and bottom of range
if ~isempty(Fnum), 
	if Fnum==0, Fnum=figure; else figure(Fnum); end
	clf; hold on;
	plot([min([min(Rhat) min(Robs)]) max([max(Rhat) max(Robs)])],[min([min(Rhat) min(Robs)]) max([max(Rhat) max(Robs)])],'r--','LineWidth',1.4);
	plot(Robs,Rhat,'ko','MarkerFaceColor',[.8 .6 .2]); 
	Rhat=raylcdf(Dlist,Bhat(1)); 
	plot(Robs,Rhat,'ko','MarkerFaceColor',[.7 .5 .1]); end

%REDUCE range below TOL [assumption: Btrue~=mean(Bhat)]
Niterate=0;
while and(diff(Bhat)>TOL,Niterate<Ncrit), 
	Niterate=Niterate+1;
	Rhat=raylcdf(Dlist,mean(Bhat));
	Bhat((sum(Rhat-Robs)<0)+1)=mean(Bhat); end
if Niterate==Ncrit, warning('Search terminated before reaching criterion'); end %#ok<WNTAG>

%PLOT result for mean(Bhat) along with results from range
if ~isempty(Fnum), 
	plot([min([min(Rhat) min(Robs)]) max([max(Rhat) max(Robs)])],[min([min(Rhat) min(Robs)]) max([max(Rhat) max(Robs)])],'k--','LineWidth',1);
	plot(Robs,Rhat,'ko','MarkerFaceColor',[.4 .4 .8]); end

%VERIFY that computed rT is not larger than Rmax
if isempty(Rmax), rT=raylinv(DesiredHitRate,mean(Bhat)); else rT=min([raylinv(DesiredHitRate,mean(Bhat)) Rmax]); end
rTrange=raylinv(DesiredHitRate,Bhat);

