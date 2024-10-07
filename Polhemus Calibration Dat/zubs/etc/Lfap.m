% Lfap.m
%
% USAGE:
% [FAPhat Lcell Lcells DWmat]=Lfap(Din,Vcell,varargin); %nan's are treated as missing data
%
% INPUTS
%      Din: data (COL1) and trial numbers or times of each datum (COL2)
%    Vcell: 3 cells: frequencies, amplitudes, and phases to compute the prob of
% Optional:
%    Trend: 0-2; indicates order (+1) of Taylor expansion describing trend
%      0:no trend [default], 1:unknown offset, 2:unknown linear trend
%      **Note trend can be given as first optional arg without naming it
%    Epoch: Length of epoch (trials or seconds); used to give experiment-sensible phase numbers. If empty, phase is given in radians
%   Fnum: Number of first of two figures for plotting (second fig is Fnum+1)
%   Fprime: Theoretical frequency (used for plotting, and to estimate Alike and Plike; peak of Flike is used if no Fprime is provided)
%   Aprime: Theoretical amplitude (used for plotting, and output is given as Amplitude gain if Aprime is provided)
%   Pprime: Theoretical phase (used for plotting)
%
% OUTPUTS
%  FAPhat: struct containing [max mean SD] of f, a, and p values (overall
%          on row1, separate sources on subsequent rows)
%   Lcell: {1}=[Lf Flist], {2}=[La Alist], {3}=[Lphi Plist]
%  Lcells: same as Lcell, but the first S columns are log-probabilities by subject (and last col is Flist, Alist, or Plist)
%   DWmat: [plotdat wm t]
% *** Note that max values are interpolated
%
% OUTPUT PLOTS 
%  F1.1: Orange.[data]        F1.2-F2.2 k.[log(p)]
%        Grey.[windowedmean]            k--[prime]
%        k--[pert]                      o[best-fit]
%        Grey--[best-fit]
%
% %%%%%%%%%%%%%%%%
% %%% EXAMPLES %%%
% Fprime=5; Aprime=2; Pprime=-pi/12;
% offset=.5*Aprime*[-1 0 1]; slope=5; SD=1.75*(2*Aprime);
% T=517; t=[0:T-1]'/T;
% Flist=linspace(1,12,12*14-(14-1));
% Alist=Aprime*linspace(-1.75,1.75,302);
% Plist=linspace(-pi,pi,64);
% Vin=cell(3,1); Vin{1}=Flist; Vin{2}=Alist; Vin{3}=Plist;
% noise=SD*(rand([length(t) 5])-.5);
%
% %COS SIGNAL
% dxn=Aprime*cos(t*2*pi*Fprime+Pprime)+noise(:,1); 
% Dnow=cell(1,2); Dnow{1}=[dxn(:,1) t];
% [FAPhat Lcell]=Lfap(Dnow,Vin,0,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
%
% %SIN SIGNAL (MULTIPLE DATASETS)
% dxn=Aprime*sin(t*2*pi*Fprime+Pprime)*[1 1 1 1 1]+noise; 
% Dnow=cell(5,2); for s=1:5, Dnow{s,2}=[dxn(:,s) t]; end
% [FAPhat Lcell]=Lfap(Dnow,Vin,0,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
% 
% %SIN+OFFSET
% dxn=Aprime*sin(t*2*pi*Fprime+Pprime)+noise(:,2); 
% Dnow=cell(1,2); Dnow{2}=[dxn+offset(1) t];
% [FAPhat Lcell]=Lfap(Dnow,Vin,1,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
%
% %SIN+OFFSET (MULTIPLE DATASETS; MISSING DATA, UNEQUALLY-SIZED DATASETS, UNEQUAL OFFSETS)
% for s=1:3, inds=randperm(T); inds=sort(inds(1:round(T/10)+randi([-10 50])))'; Dnow{s,2}=[dxn(inds)+offset(s) t(inds)]; end
% [FAPhat Lcell]=Lfap(Dnow,Vin,0,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
%
% %SIN+SLOPE (uses 'Epoch' input)
% dxn=slope*t+Aprime*sin(t*2*pi*Fprime+Pprime)+noise(:,2); 
% Dnow=cell(1,2); Dnow{2}=[dxn t];
% [FAPhat Lcell]=Lfap(Dnow,Vin,1,'Fnum',1,'Epoch',T,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
% [FAPhat Lcell]=Lfap(Dnow,Vin,2,'Fnum',1,'Epoch',T,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
%
% %SIN+OFFSET+SLOPE
% dxn=offset+slope*t+Aprime*sin(t*2*pi*Fprime+Pprime)+noise(:,3); 
% Dnow=cell(1,2); Dnow{2}=[dxn t];
% [FAPhat Lcell]=Lfap(Dnow,Vin,2,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
%
% %SIN+OFFSET+SLOPE (MISSING DATA)
% inds=randi(T,[round(4*T/5) 1]); Dnow{2}(inds,1)=nan;  
% [FAPhat Lcell]=Lfap(Dnow,Vin,2,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
%
% %SIN+OFFSET+SLOPE+QUADRATIC
% dxn=offset+slope*t+25*(t-.3).^2+Aprime*sin(t*2*pi*Fprime+Pprime)+noise(:,3); 
% Dnow=cell(1,2); Dnow{2}=[dxn t];
% [FAPhat Lcell]=Lfap(Dnow,Vin,0,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
% [FAPhat Lcell]=Lfap(Dnow,Vin,1,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
% [FAPhat Lcell]=Lfap(Dnow,Vin,2,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
% [FAPhat Lcell]=Lfap(Dnow,Vin,3,'Fnum',1,'Aprime',Aprime,'Fprime',Fprime,'Pprime',Pprime);
%
% teh wrote it (10.20.11)
% updated (4.11.12) to include expanded commenting, optional Epoch variable, use of logsum function, and updated plotting
% updated (4.16.12) to include compensation for quadratic trend, eliminate possibility of negative maxAmp, and made changes to speed primary loop
% updated (4.24.12) to fix issue where  maxAmp is found at zero, as interpolation between two identical peaks at +- AmpMax; added plots of best-fit signal overlaid on data (and reference signal with negative amplitude when Sprime values are input
% updated (4.27.12) to fix plotting issue when Fprime is not provided; include plotting of average binned data if desired, added examples with missing data
% updated (3 14.15) to replace slow shortestCI.m with faster CIp.m for finding confidence intervals
% updated (8.07.15) updated platted to include flags to add FIT, signal, and WMplot plots to the dataplot


function [FAPhat Lcell Lcells DWmat]=Lfap(Din,Vcell,varargin) %#ok<*TRYNC>
if iscell(Din), Dcell=Din; S=size(Dcell,1);
	if size(Dcell,2)==1, Dcell(:,2)=cell(S,1); end
else S=size(Din,3); Dcell=cell(S,2);
	if and(~(size(Dcell,4)>1),~(size(Dcell,2)>2)),
		quans=questdlg('sin or cos data','','sin','cos','sin'); end
	for s=1:S,
		if size(Dcell,4)>1,
			if all(isnan(Din(:,1,s,1))), Dcell{s,1}=[];
			else Dcell{s,1}=Din(:,:,s,1); end
			Dcell{s,2}=Din(:,:,s,2);
		elseif size(Dcell,2)>2,
			if all(isnan(Din(:,1,s))), Din{s,1}=[];
			else Dcell{s,1}=Din(:,1:2,s); end
			Dcell{s,2}=Din(:,3:4,s);
		else if strcmp(quans,'sin'), Dcell{s,2}=Din(:,:,s);
			else Dcell{s,1}=Din(:,:,s); end, end, end, end; clear Din
base='common';
trend=0;
if nargin==1, Flist=linspace(0,15,301)'; Alist=linspace(-1.25,1.25,301); Plist=linspace(-pi,pi,161)';
elseif and(nargin==2,isscalar(Vcell)), trend=Vcell; Flist=linspace(0,15,301)'; Alist=linspace(-1.25,1.25,301); Plist=linspace(-pi,pi,161)'; clear Vcell
else if length(Vcell{1})==2, Flist=linspace(0,Vcell{1}(1),Vcell{1}(2))';
	elseif xor(size(Vcell{1},1)==1,size(Vcell{1},2)==1), Flist=Vcell{1}(:); else Flist=Vcell{1}; end
	if length((Vcell{2}))==1, Alist=linspace(-1.75,1.75,Vcell{2});
	else if xor(size(Vcell{2},1)==1,size(Vcell{2},2)==1), Alist=ones(S,1)*Vcell{2}(:)';
		else if size(Vcell{2},1)==S, Alist=Vcell{2}; elseif size(Vcell{2},2)==S, Alist=Vcell{2}'; else error('Vin{2} must be a scalar, vector, or an SxM matrix'); end, end, end
	if length(Vcell{3})==1, Plist=linspace(-pi,pi,Vcell{3})'; else Plist=Vcell{3}; end, end
if and(ismatrix(Flist),size(Flist,1)==1), Flist=Flist(:); end
if and(ismatrix(Alist),size(Alist,2)==1), Alist=Alist(:)'; end
if and(ismatrix(Plist),size(Plist,1)==1), Plist=Plist(:); end
Lcells=cell(3,1);


%Variable (optional) input arguments
MASH=1; FIT=[]; SIG=[]; FORCE=0; WMplot=0;
Fnum=[]; Fprime=nan; Aprime=nan; Pprime=nan; Epoch=nan; Tmin=[]; Tmax=[]; CI=[];
if rem(length(varargin),2)==0, va=varargin;	
elseif and(isnumeric(varargin{1}),all(size(varargin{1}))==1), va=varargin(2:end); trend=varargin{1}; 
else error('First variable argumenst must be the trend number or a variable descriptor string'); end
for n=1:length(va)/2,
	switch lower(va{1}),
		case 'force', FORCE=1;
		case 'fit', FIT=va{2};
		case 'signal', SIG=-va{2};
		case 'wmplot', WMplot=va{2};
		case 'trend', if and(isnumeric(va{2}),all(size(va{2}))==1), trend=va{2}; if or(trend<0,trend>3), error('Trend must be between 0 and 3'); end, else error('Incorrect assignment to property ''trend'''); end
		case 'epoch', if and(isnumeric(va{2}),all(size(va{2}))==1), Epoch=va{2}; else error('Incorrect assignment to property ''Epoch'''); end
		case 'fignum', if or(isempty(va{2}),and(isnumeric(va{2}),all(size(va{2}))==1)), Fnum=va{2}; else error('Incorrect assignment to property ''Fnum'''); end
		case 'fnum', if or(isempty(va{2}),and(isnumeric(va{2}),all(size(va{2}))==1)), Fnum=va{2}; else error('Incorrect assignment to property ''Fnum'''); end
		case 'ci', if or(isempty(va{2}),and(isnumeric(va{2}),all(size(va{2}))==1)), CI=va{2}; else error('Incorrect assignment to property ''Fnum'''); end
		case 'mash', if and(isnumeric(va{2}),all(size(va{2}))==1), MASH=va{2}; else error('Incorrect assignment to property ''MASH'''); end
		case 'fprime', if and(isnumeric(va{2}),all(size(va{2}))==1), Fprime=va{2}; else error('Incorrect assignment to property ''Fprime'''); end
		case 'aprime', if and(isnumeric(va{2}),all(size(va{2}))==1), Aprime=va{2}; else error('Incorrect assignment to property ''Aprime'''); end
		case 'pprime', if isnumeric(va{2}), if max(size(va{2}))==1, Pprime=va{2}; end, else error('Incorrect assignment to property ''Pprime'''); end
		case 'phiprime', if isnumeric(va{2}), if max(size(va{2}))==1, Pprime=va{2}; end, else error('Incorrect assignment to property ''PHIprime'''); end
	end, if length(va)>3, va=va(3:end); end, end
if and(size(Alist,2)==1,S>1), Alist=ones(S,1)*Alist; end
if and(length(Aprime)==1,S>1), Aprime=Aprime*ones(S,1); end
if isempty(CI), CI=.95; end
if SIG==0, SIG=1; end
	
%Initialize output variables
FAPhat=struct('key',['peak / mean / SD / CImin / CImax'],'f',[],'amp',struct('real',[],'gain',[]),'phase',struct('rad',[],'trial',[]));
Prad=Plist; if ~isnan(Pprime), Pprad=Pprime; end
%if and(~isnan(Epoch),~isnan(Fprime)), Plist=Epoch/Fprime*Plist/(2*pi); 
%	if ~isnan(Pprime), Pprime=Epoch/Fprime*Pprime/(2*pi); end, end
Lout=zeros(length(Flist),size(Alist,2),length(Plist)); Ltmp=Lout;
for s=1:S, ALnow=Alist(s,:);
	if isempty(Dcell{s,2}), D=Dcell{s,1}(:,1)'; fnind=1; else D=Dcell{s,2}(:,1)'; fnind=2; end
	T=length(D); if size(Dcell{s,fnind},2)==1, t=linspace(0,1,T+1)'; t=t(1:end-1); else t=Dcell{s,fnind}(:,2)'; end
	indnnan=find(~isnan(D)); D=D(indnnan); t=t(indnnan); T=length(indnnan);
	omeg=(2*pi*Flist)*t;
	Vd=zeros(trend+1,1); Vd(1)=mean(D.^2);
	if trend>0, mD=mean(D); Vd(2)=-mD^2;
		if trend>1, mT=mean(t); mT2=mean(t.^2); mDT=mean(D.*t);
			Vt=mT2-mT^2; Vdt=(mDT-mD*mT); Vd(3)=-Vdt^2/Vt;
			if trend==3, mDT2=mean(D.*(t.^2));
				Vtu=mean(t.^3)-mT2*mT; Vt2=Vtu/Vt; 
				VT2=(mean(t.^4)-mT2^2)-Vtu^2/Vt;
				Vdt2=(mDT2-mD*mT2)-(mDT-mD*mT)*Vt2; Vd(4)=-Vdt2^2/VT2; end, end, end
	tmp1=zeros(length(Flist),trend+1);
	tmp2=zeros(length(Flist),trend+1);
	for pind=1:length(Plist), Pnow=Prad(pind);
		if fnind==1, C=cos(omeg+Pnow); else C=sin(omeg+Pnow); end
		tmp1(:,1)=mean(C.^2,2);
		tmp2(:,1)=mean(C.*(ones(size(Flist))*D),2);
		if trend>0, mC=mean(C,2);
			tmp1(:,2)=-mC.^2; 
			tmp2(:,2)=-mD*mC;
			if trend>1, mCT=mean(C.*(ones(size(Flist))*t),2); Vct=mCT-mT*mC;
				tmp1(:,3)=-Vct.^2/Vt; 
				tmp2(:,3)=-Vct*Vdt/Vt;
				if trend==3, mCT2=mean(C.*(ones(size(Flist))*(t.^2)),2);
					Vct2=(mCT2-mC*mT2)-(mCT-mC*mT)*Vt2;
					tmp1(:,4)=-Vct2.^2/VT2;
					tmp2(:,4)=-Vct2*Vdt2/VT2; end, end, end
		Ltmp(:,:,pind)=-(T-trend)*...
			log10(sum(Vd)+sum(tmp1,2)*(ALnow.^2)-2*sum(tmp2,2)*ALnow)/2; end
	
	inclist=[diff(Flist(1:2)) diff(ALnow(1:2)) diff(Plist(1:2))];
	Lout=Lout+Ltmp;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Best-fit parameters %%%
	if S>1, [FAPhat,~,~,~,Lcell]=getFAP(FAPhat,s+1,Flist,ALnow,[Plist(:) Prad(:)],[Fprime Aprime(s) Pprime],Epoch,Ltmp,base,inclist,CI); 
		for rownow=1:3, Lcells{rownow}(:,s:s+1)=Lcell{rownow}; end, end, end
[FAPhat Alist Plist Pprime Lcell]=getFAP(FAPhat,1,Flist,ALnow,[Plist(:) Prad(:)],[Fprime Aprime(s) Pprime],Epoch,Lout,base,inclist,CI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compile avg binned data %%%
if or (~isempty(Fnum),nargout>2), %set DWmat output (also used for plotting)
	Nbins=[0 0];
	if or(S==1,~MASH), Dcell=Dcell(end,:); end
	for fnum=1:2, bins=0; for s=1:size(Dcell,1),
			if size(Dcell{s,fnum},2)<2, t=[0 1]'; else t=Dcell{s,fnum}(:,2); end
			
			t
			
			T=size(Dcell{s,fnum},1); bins=bins+T; Tmin=min([Tmin; t(:)]); Tmax=max([Tmax; t(:)]); end
		
		
		if bins>=500, bins=101; else bins=round(bins/(1+sum(.25*[bins>100 bins>200 bins>300 bins>350 bins>400]))); end
		Nbins(fnum)=bins; end
	DWmat=zeros(max(Nbins),3,2);
	for fnum=1:2, clear Bcell Bmat
		if ~isempty(Dcell{1,fnum}),
		for s=1:size(Dcell,1), t=[]; D=[];
			if size(Dcell{s,fnum},2)==1, 1
				t=linspace(0,1,size(D,1)+1)'; t=t(1:end-1); D=Dcell{s,fnum}/Aprime(s)
			elseif 2
				size(Dcell{s,fnum},2)==2, t=Dcell{s,fnum}(:,2); D=Dcell{s,fnum}(:,1)/Aprime(s)
	end
			if s==1, [Bcell Bmat]=binner(D,t,[Tmin Tmax max(Nbins)]); 
			else [Bcell]=binner(D,t,Bcell,Bmat); end, end
		for bt=1:size(Bmat,1), Bmat(bt,:)=[nanmean(Bcell{bt}) mean(Bmat(bt,:))]; end %note: nanmean([])=nan
		
		if ~isnan(Fprime), wm=windowmean(Bmat(:,1),round(length(Bmat(:,1))/(Fprime*range(t))/2),[]); %errors may arise when Flike is weakly peaked and Fprime is not given
		else wm=windowmean(Bmat(:,1),round(length(Bmat(:,1))/(FAPhat.f(1)*range(t))/2),[]); end
		DWmat(:,:,fnum)=[mean(Aprime)*[Bmat(:,1) wm(:)] Bmat(:,2)]; end, end, end

% figure; clf; meshc(Plist,Alist,real(Ltmp));
% ylabel('Amp','FontName','Times','FontWeight','Bold','FontSize',14);
% xlabel('\phi','FontWeight','Bold','FontSize',16); r=axis;
% axis([Plist(1) Plist(end) Alist(1) Alist(end) r(5:6)]); colormap(hot)

%%%%%%%%%%%%%%%%%%%%%
%%% LFAP PLOTTING %%%
if ~isempty(Fnum), %Plotting ONLY -- from here to end of Lfap function
	Dcol=[.2 .6 .8;.8 .6 .2];
	% OPTIONAL PLOT
		if Fnum==0, Fnum=figure; else figure(Fnum); end
	clf; subplot(2,1,1); hold on
	
	DWmat
	
	%Plot DATA
	for fnum=1:2, if fnind==fnum,
			if S==1, 1
				plot(DWmat(:,3,fnum),DWmat(:,1,fnum),'ko','Color',Dcol(fnum,:),'MarkerSize',8);
			else 2
				plot(DWmat(:,3,fnum),DWmat(:,1,fnum),'ko','MarkerFaceColor',Dcol(fnum,:),'MarkerSize',6); end
			if WMplot==1,
				try	plot(DWmat(:,3,fnum),DWmat(:,2,fnum),'-','Color',.65*[1 1 1],'LineWidth',3); end, end, end, end
	
	if any(fnind)==1, %sin(fnind=2) or cos(=1)-phase inducer/data
		if and(~isnan(Fprime),and(~isnan(Aprime),~isnan(Pprime))),
			%SIG indicates whether to plot the (neg) pert signal
			if ~isempty(SIG), plot(t,SIG*Aprime(end)*cos(2*pi*Fprime*t+Pprad),'k--','LineWidth',1.1); end, end
		if ~isempty(FIT), %FIT indicates whether to plot the fitted sinusoid
			if or(~FORCE,isnan(Fprime)), plot(t,FIT*FAPhat.amp.real(1)*cos(2*pi*FAPhat.f(1)*t+FAPhat.phase.rad(1)),'--','Color',.45*[1 1 1],'LineWidth',1.1);
			elseif FORCE, plot(t,FIT*FAPhat.amp.real(1)*cos(2*pi*Fprime*t+FAPhat.phase.rad(1)),'--','Color',.45*[1 1 1],'LineWidth',1.1); end, end, end
	if any(fnind)==2,
		if and(~isnan(Fprime),and(~isnan(Aprime),~isnan(Pprime))),
			if ~isempty(SIG), plot(t,SIG*Aprime(end)*sin(2*pi*Fprime*t+Pprad),'k--','LineWidth',1.1); end, end
		if ~isempty(FIT),
			if or(~FORCE,isnan(Fprime)), plot(t,FIT*FAPhat.amp.real(1)*sin(2*pi*FAPhat.f(1)*t+FAPhat.phase.rad(1)),'--','Color',.45*[1 1 1],'LineWidth',1.1);
			elseif FORCE, plot(t,FIT*FAPhat.amp.real(1)*sin(2*pi*Fprime*t+FAPhat.phase.rad(1)),'--','Color',.45*[1 1 1],'LineWidth',1.1); end, end, end
	r=axis; plot(r(1:2),[0 0],'k-','LineWidth',.7)
	xlabel(['\Deltat'],'FontName','Arial','FontSize',16)
	ylabel('\Deltay','FontName','Arial','FontSize',16)
	
	%Plot log-Fdist
	subplot(2,1,2); hold on
	plot(Flist,Lcell{1}(:,1),'k.-'); box off; r=axis;
	if ~isnan(Fprime), plot(Fprime*[1 1],[r(3) 0],'k--'); end
	plot(FAPhat.f(1),r(3),'ko','MarkerSize',5,'MarkerFaceColor',[.7 .7 1]);
	try axis([Flist(1) Flist(end) r(3) .025*abs(r(3))]); end
	xlabel(['frequency'],'FontName','Arial','FontSize',16);
	ylabel('\Delta\Lambda','FontSize',14);
	
	%Plot log-Adist
	figure(Fnum+1); clf
	subplot(2,1,1); hold on
	plot(Alist,Lcell{2}(:,1),'k-','LineWidth',2); box off; r=axis;
	if isnan(Aprime), Alabel='amplitude'; anow=FAPhat.amp.real(1);
	else plot([1 1],[r(3) 0],'k--'); Alabel='amplitude gain'; anow=FAPhat.amp.gain(1); end
	plot(Alist(Alist<=0),Lcell{2}(Alist<=0,1),'.','Color',.7*[1 1 1]);
	plot(Alist(Alist>=0),Lcell{2}(Alist>=0,1),'k.');
	plot(anow,r(3),'ko','MarkerSize',5,'MarkerFaceColor',[.7 .7 1]);
	plot([0 0],r(3:4),'k:'); try axis([Alist(1) Alist(end)  r(3) .025*abs(r(3))]); end
	xlabel(Alabel,'FontName','Arial','FontSize',16)
	ylabel('\Delta\Lambda','FontSize',14);
	
	%Plot log-Pdist
	subplot(2,1,2); hold on
	if isnan(Epoch), Plabel='phase (rad)'; else Plabel='phase (t)'; end
	plot(Plist,Lcell{3}(:,1),'k-','LineWidth',2);
	if length(Plist)<250, LW=4; MY='-'; else LW=1; MY='.'; end
	ind0=find(Plist>=0,1); ind=[find(isnear(Lcell{3}(1:ind0-1,1),min(Lcell{3}(1:ind0-1,1))),1) ind0+find(isnear(Lcell{3}(ind0:end,1),min(Lcell{3}(ind0:end,1))),1)-1];	
	plot(Plist(1:ind(1)+1),Lcell{3}(1:ind(1)+1,1),MY,'Color',.7*[1 1 1],'LineWidth',LW);
	plot(Plist(ind(2)-1:end),Lcell{3}(ind(2)-1:end,1),MY,'Color',.7*[1 1 1],'LineWidth',LW);
	plot(Plist(ind(1)+1:ind(2)-1),Lcell{3}(ind(1)+1:ind(2)-1,1),MY,'Color','k','LineWidth',LW); box off; r=axis;
	if ~isnan(Pprime), if Pprime~=0, plot(Pprime*[1 1],[r(3) 0],'k--'); end, end
	if ~isnan(Epoch), plot(FAPhat.phase.trial(1,1),r(3),'ko','MarkerSize',5,'MarkerFaceColor',[.7 .7 1]);
	else plot(FAPhat.phase.rad(1,1),r(3),'ko','MarkerSize',5,'MarkerFaceColor',[.7 .7 1]); end
	plot([0 0],r(3:4),'k:'); try axis([Plist(1) Plist(end)  r(3) .025*abs(r(3))]); end
	if isnan(Epoch), set(gca,'XTick',-pi:pi/2:pi);
		set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'}); end
	xlabel(Plabel,'FontName','Arial','FontSize',16)
	ylabel('\Delta\Lambda','FontSize',14); end, end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  ZUBS  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%
%%% getFAP.m %%%
function [FAP Alist Plist Pprime Lcell]=getFAP(FAP,rowNum,Flist,Alist,Plist,Primes,Epoch,Lnow,base,inclist,CI)
Fprime=Primes(1); Aprime=Primes(2); Pprime=Primes(3); Prad=Plist(:,2); Plist=Plist(:,1); 

%Frequency distribution
Lf=logsum(Lnow,[2 3],inclist(2:3),'baseIn',base); 
Lf=Lf(:)-max(Lf);
FAP.f(rowNum,1)=logpeakval([Lf(:) Flist(:)]); 
%[Lf(:) Flist(:)]
try [V M]=vardef(Lf,Flist(:)); FAP.f(rowNum,2:3)=[M sqrt(V)]; end
try FAP.f(rowNum,4:5)=CIp(Flist(:),Lf,CI); end

%Shrink Lnow based on Fprime (if unavailable, Fest)
if ~isnan(Fprime), ind=find(abs(Fprime-Flist)==min(abs(Fprime-Flist)),1); else ind=find(abs(FAP.f(rowNum,1)-Flist)==min(abs(FAP.f(rowNum,1)-Flist)),1); end
Ltmp=squeeze(Lnow(ind,:,:));

%Amplitude distribution
La=logsum(Ltmp,2,inclist(3),'baseIn',base); La=La(:)-max(La(:));
FAP.amp.real(rowNum,1)=logpeakval([La(Alist>=0) Alist(Alist>=0)']); try [V M]=vardef(La(Alist>=0),str8n(Alist(Alist>=0))); FAP.amp.real(rowNum,2:3)=[M sqrt(V)]; end
try FAP.amp.real(rowNum,4:5)=CIp(str8n(Alist(Alist>=0)),La(Alist>=0),CI); end
if ~isnan(Aprime), Alist=Alist/Aprime;
	FAP.amp.gain(rowNum,1)=logpeakval([La(Alist>=0) Alist(Alist>=0)']);
	try [V M]=vardef(La(Alist>=0),str8n(Alist(Alist>=0))); FAP.amp.gain(rowNum,2:3)=[M sqrt(V)]; end
	try FAP.amp.gain(rowNum,4:5)=CIp(str8n(Alist(Alist>=0)),La(Alist>=0),CI); end, end

%Phase distribution
Lp=logsum(Ltmp,1,inclist(2),'baseIn',base); Lp=Lp(:)-max(Lp(:));
ind0=[find(Prad>=-pi/2,1) find(Prad>=pi/2,1)]; if length(ind0)==1, ind=find(Lp==max(Lp)); else ind=ind0(1)+find(Lp(ind0(1):ind0(2))==max(Lp(ind0(1):ind0(2))))-1; end
FAP.phase.rad(rowNum,1)=logpeakval([Lp(:) Prad(:)],ind);
ind0=find(Plist>=0,1); 
ind=[find(Lp(1:ind0-1)==min(Lp(1:ind0-1)),1) ind0+find(Lp(ind0:end)==min(Lp(ind0:end)),1)-1];
%[str8n(Plist(ind(1):ind(2))) str8n(Prad(ind(1):ind(2)))]
try [V M]=vardef(Lp(ind(1):ind(2)),str8n(Prad(ind(1):ind(2)))); %vardef(Lp(ind(1)+1:ind(2)-1),str8n(Prad(ind(1)+1:ind(2)-1)));

FAP.phase.rad(rowNum,2:3)=[M sqrt(V)]; end
try FAP.phase.rad(rowNum,4:5)=CIp(str8n(Prad(ind(1)+1:ind(2)-1)),Lp(ind(1)+1:ind(2)-1),CI); end
%Convert Plist to experimentally sensible units
if ~isnan(Epoch),
	if ~isnan(Fprime), Plist=Epoch/Fprime*Prad/(2*pi);
		if ~isnan(Pprime), Pprime=Epoch/Fprime*Pprime/(2*pi); end
	else Plist=Epoch/FAP.f(rowNum,1)*Prad/(2*pi);
		if ~isnan(Pprime), Pprime=Epoch/FAP.f(rowNum,1)*Pprime/(2*pi); end, end
	try FAP.phase.trial(rowNum,1)=logpeakval([Lp(ind(1)+1:ind(2)) Plist(ind(1)+1:ind(2))]); end
	
% 	[Plist(ind(1):ind(2)) Prad(ind(1):ind(2))]
% 	Fprime
% 	Aprime
% 	Pprime
% 	Epoch
	
	try [V M]=vardef(Lp(ind(1):ind(2)),str8n(Plist(ind(1):ind(2)))); %[V M]=vardef(Lp(ind(1)+1:ind(2)-1),str8n(Plist(ind(1)+1:ind(2)-1))); 
	FAP.phase.trial(rowNum,2:3)=[M sqrt(V)]; end
	try FAP.phase.trial(rowNum,4:5)=CIp(str8n(Plist(ind(1)+1:ind(2)-1)),Lp(ind(1)+1:ind(2)-1),CI); end, end

Lcell=cell(4,1);
Lcell{1}=[Lf(:) Flist(:)]; Lcell{2}=[La(:) Alist(:)]; Lcell{3}=[Lp(:) Plist(:)]; Lcell{4}=Lnow; end


% logpeakval.m
%
% usage: [logPout]=logpeakval(logPvec);
%
% Computes interpolated peak location given index of sampled peak
%
% INPUTS
% logPvec: [likelihoods abscissaValues]
% optional:
%     ind: index or indices of peak
%     epn: exponent for weighting (bet .25 and .5 works best)
%
% OUTPUT
% peakval: abscissa value corresponding to interpolated peak
%
% teh wrote it. [10.12.11] ** DO NOT DISTRIBUTE

function [peakval]=logpeakval(logPvec,ind,epn)
if all(rem(logPvec(:,1),1)==0), logPvec(:,1)=log(logPvec(:,1)/sum(logPvec(:,1))); end
if any(size(logPvec)==1), error('logPvec input must be a 2-column matrix: C1=p; C2=x'); end
if nargin<3, epn=.375; ind=find(logPvec(:,1)==max(logPvec(:,1))); end
if nargin==2, if all(rem(ind,1)==0), epn=.375; else epn=ind; ind=find(logPvec(:,1)==max(logPvec(:,1))); end, end
if length(ind)==1,
	if or(ind==1,ind==length(logPvec)), peakval=logPvec(ind,2);
	else logPin=logPvec(ind-1:ind+1,1);
		d1=logPin(2)-logPin(1);
		d2=logPin(2)-logPin(3);
		peakval=(d1^epn*logPvec(ind+1,2)+d2^epn*logPvec(ind-1,2))/(d1^epn+d2^epn); end
else peakval=mean(logPvec(ind,2)); end, end


% logsum.m
%
% usage: [Lout]=logsum(Lin,dim);
%        [Lout]=logsum(Lin,dim,trapinc);
%        [Lout]=logsum(Lin,dim,'trapinc',trapinc,'baseIn',base);
%
% Computes binary log of a sum of the log2 of matrix/vector elements.
% Used primarily when summing over one dimension of a multidimensional
% log-probability distribution
% Based on the identity:
%    log(sum(a))=log(a1)+log(1+sum(10.^(log(a(2:end))-log(a1))))
%
% INPUTS
%  Lin: multidimensional input distribution of the form:
%           Lin=log10(Pmat); [must be log10]
%     dim: dimension that is to be eliminated from Pmat
%           Note: if dim is a vector, all dimensions named in the vector are
%           eliminated. The order is irrelevant.
% OPTIONAL
% trapinc: Used if the sum is to be done via trapezoidal integration
%           Gives the increment between axis values [DEFAULT=1]
%           Note: size(trapinc), if trapinc is given, must be same as size(dim).
%  baseIn: gives the base used to create logPmat; 'binary, 'common', or 'natural'
%          DEFAULT='natural' [Note legacy default='common']
%
% OUTPUT
% Lout: Lout=log(sum(exp(logPmat)));
%
% EXAMPLE: test values
% lp=log2([.1 .2 .3 .4 .4 .3 .2 .1]);
% logsum(lp,2,'baseIn','binary')
%
% EXAMPLE: test dims
% lp=log(rand([4 2 6 1 3 1 5]));
% size(logsum(lp,[1 3 4],'baseIn','natural'))
%
% EXAMPLE: test values and dims
% mu=[2 -2]; S=[3 0;0 2]; x=linspace(-10,10,41); y=linspace(-10,10,41);
% [x2 x1]=meshgrid(y,x);
% X=[x1(:) x2(:)];
% P=reshape(mvnpdf(X,mu,S),length(x),length(y)); lP=log(P);
% figure; meshc(x1,x2,P); colormap(bone)
% xlabel('X','FontSize',16,'FontWeight','Bold')
% ylabel('Y','FontSize',16,'FontWeight','Bold')
% zlabel('prob density','FontSize',16,'FontWeight','Bold')
% figure; subplot(1,2,1); hold on
% plot(x,logsum(lP,2,x(2)-x(1)),'k.'); xlabel('X-axis','FontSize',16,'FontWeight','Bold')
% plot(x,log(mvnpdf(x',mu(1),S(1,1))),'--','LineWidth',1.1)
% subplot(1,2,2); hold on
% plot(y,logsum(lP,1,y(2)-y(1)),'k.'); xlabel('Y-axis','FontSize',16,'FontWeight','Bold')
% plot(y,log(mvnpdf(y',mu(2),S(2,2))),'--','LineWidth',1.1)
%
% EXAMPLE:
% mu=[1 2 3]; S=[3 2 1]; x=linspace(-10,10,100);
% [x2 x1 x3]=meshgrid(x',x',x');
% X=[x1(:) x2(:) x3(:)];
% P=reshape(mvnpdf(X,mu,S),100,100,100); lP=log10(P);
% figure; hold on
% plot(x,logsum(lP,[1 3],x(2)-x(1)*[1 1],'baseIn','common'),'k.')
% plot(x,log10(trapz(trapz(P,3)*(x(2)-x(1)),1)*(x(2)-x(1))),'--','LineWidth',1.1)
% xlabel('Y-axis','FontSize',16,'FontWeight','Bold')
%
% teh wrote it. [10.12.11] ** DO NOT DISTRIBUTE

function Lout=logsum(Lin,dim,varargin)
if all(size(Lin)==1), Lout=Lin;
else
	%%%%%%%%%%%%%%%%%%%%%%
	%%% Set the inputs %%%
	baseIn='natural'; trapinc=1; %defaults
	if nargin==1,
		if ~and(ndims(Lin)==2,min(size(Lin))==1), %#ok<ISMAT>
			error('Lin must be a vector when it is the only input to logsum');
		else dim=1; end
	else %nargin>1
		if ischar(dim), baseIn=dim; dim=1; end
		trapinc=ones(size(dim));
		
		if ~isempty(varargin), va=varargin;
			if rem(length(va),2)==1,
				if ischar(va{1}), baseIn=va{1}; else trapinc=va{1}; end
				va=va(2:end); end
			if length(va)==2, %disambiguates inputs when no namestrings are given (dim must precede trapinc)
				if ischar(va{1}), if or(or(strcmp(lower(va{1}),'binary'),strcmp(lower(va{1}),'natural')),strcmp(lower(va{1}),'common')), baseIn=va{1}; trapinc=va{2}; va=[]; end
				elseif ischar(va{2}), if or(or(strcmp(lower(va{2}),'binary'),strcmp(lower(va{2}),'natural')),strcmp(lower(va{2}),'common')), trapinc=va{1}; baseIn=va{2}; va=[]; end
				elseif and(~ischar(va{1}),~ischar(va{2})), dim=va{1}; trapinc=va{2}; va=[]; end, end
			Nva=length(va)/2+1;
			while length(va)>1, Nva=Nva-1;
				if Nva==0, error('unmatched input string'); end
				switch lower(va{end-1}),
					case 'trapinc', trapinc=va{end}; va=va(1:end-2);
					case 'seplist', trapinc=va{end}; va=va(1:end-2);
					case 'basein', baseIn=va{end}; va=va(1:end-2);
					case 'base', baseIn=va{end}; va=va(1:end-2); end, end, end, end
	if and(ndims(Lin)==2,min(size(Lin))==1), Lin=Lin(:); end %#ok<ISMAT>
	if any(dim>ndims(Lin)), error('dim input exceeds dimensionality of Lin'); end
	if length(dim)~=length(trapinc), error('must include trapinc spacings for all integrated dimensions (or none)'); end
	
	%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Change to BINARY %%%
	switch baseIn,
		case 'binary', Lout=Lin;
		case 'natural', Lout=Lin*log2(exp(1));
		case 'common', Lout=Lin*log2(10); end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Check for singleton dimensions %%%
	itmp=find(size(Lin)==1);
	if ~isempty(itmp), itmpdim=[];
		for m=1:length(itmp), if any(dim==itmp(m)), itmpdim(end+1)=itmp(m); end, end %#ok<AGROW>
		if ~isempty(itmpdim), itmpdim=sort(itmpdim,'ascend');
			for m=length(itmpdim):-1:1,
				trapinc=RFL(trapinc',find(dim==itmpdim(m)))';
				dim=RFL(dim',find(dim==itmpdim(m)))'; end
			Lout=rmdim(Lout,itmpdim);
			for m=1:length(itmpdim), dim(dim>itmpdim(m))=dim(dim>itmpdim(m))-1; end, end, end
	[dim I]=sort(dim(:),'ascend');
	if ~isempty(trapinc), trapinc=trapinc(I); end
		
	%%%%%%%%%%%%%%%%%
	%%% MAIN LOOP %%%
	for D=length(dim):-1:1, dnow=dim(D); nd=ndims(Lout);
		Lout=Lout+log2(trapinc(D));
		if nd>2, 
			if dnow~=nd, dimlist=1:ndims(Lout);
				Lout=permute(Lout,[dimlist(dimlist~=dnow) dnow]); end
			Lout=sort(Lout,nd,'descend');
			g=['Lout(:,:,']; gr=['1 1 ']; nm1=size(Lout,nd)-1; %#ok<NASGU>
			for n=3:nd-1, g=[g ':,']; gr=[gr '1 ']; end %#ok<AGROW>
			lp1=eval([g '1)']); lp2=eval([g '2:end)']); lprep=eval(['repmat(lp1,[' gr ' nm1])']);
		elseif nd==2,
			if dnow~=2, Lout=Lout'; end
			Lout=sort(Lout,2,'descend');
			lp1=Lout(:,1); lp2=Lout(:,2:end); lprep=Lout(:,1)*ones(1,size(Lout,2)-1);
		else nd=1; Lout=sort(Lout,'descend');
			if size(Lout,2)>1, Lout=Lout'; end
			lp1=Lout(1); lp2=Lout(2:end); lprep=lp1*ones(size(lp2)); end
		Lout=lp1+log2(1+nansum(2.^(lp2-lprep),nd)); end
	if and(nd==2,size(Lin,2)==2), Lout=Lout'; end %added specifically to deal with the case of an nx2 input matrix that is summed over rows
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Change output to baseIn %%%
	switch baseIn,
		case 'natural', Lout=Lout*log(2);
		case 'common', Lout=Lout*log10(2); end, end, end


% CIp.m
%
% Calculate p-level confidence interval of a (discretized) probability distribution
%   NOTE - this is the '2-sided' confidence interval containing p (100p%) of the
%          total probability mass
%        - This will only work for NON-UNIFORM, SINGLY-PEAKED distributions
%
% [valrange indrange interprange interpind interpvalsps]=CIp(vals,ps,p,FigNum);
% vals - abscissa values for prob dist
% ps - vector of probabilities
% p - probability mass contained in the confidence interval
% FigNum - activates plotting (does NOT use clf)
%          OPTIONAL - default   :[]
%                     'next fig':0
%
% TE Hudson
% 3.14.15

function [valrange patborders]=CIp(vals,ps,p,Ninterp)
if nargin<5, FigNum=[]; end
if nargin<4, Ninterp=[]; end
if nargin<3, p=.95; end
if nargin==2, if and(all(size(ps)==1),and(size(vals,1)>2,size(vals,2)==2)),
		p=ps; ps=vals(:,1); vals=vals(:,2); end, end
if nargin<2, error('2 inputs required'); end
if or(all(vals==0),all(ps==0)), error(['''vals'' or ''ps'' input is all zeros']); end
if p>=1, p=p/100; end
if isempty(Ninterp), Ninterp=200; end

%remove non-values and set pprime and valsprime
inds=find(and(~isnan(ps),~and(isinf(ps),sign(ps)==1)));
valsprime=vals(inds); psprime=ps(inds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forward computation %%%
if any(ps<0), cps(1)=exp(psprime(1));
	for n=2:length(psprime),
		cps(n)=exp(logsum(str8n(psprime(1:n)),1,'baseIn','natural')); end
else cps=cumsum(psprime); end

if any(diff(valsprime,1)~=diff(valsprime(1:2))),
	vals=linspace(valsprime(1),valsprime(end),length(vals));
	cpprime=interp1(valsprime,cps,vals,'pchip'); end
cpprime=cpprime/cpprime(end); fcps=cpprime;
R1range=[vals(1) vals(find(cpprime>1-p,1))]; R2range=[vals(find(cpprime<p,1,'last')) vals(end)];
Ninterp1=max([Ninterp find(cpprime>1-p,1) length(cpprime)-find(cpprime<p,1,'last')]);
vals1=linspace(R1range(1),R1range(2),Ninterp1); cps1=interp1(vals,cpprime,vals1);
vals2=linspace(R2range(1),R2range(2),Ninterp1); cps2=interp1(vals,cpprime,vals2);
dcps=str8n(cps2)*ones(1,Ninterp1)-ones(Ninterp1,1)*cps1;
dvals=diff(vals(1:2));

rangedata=nan(Ninterp1,3);
for i=1:Ninterp1, ind=find(dcps(:,i)>=p,1);
	if ~isempty(ind),
		tmpval=linspace(vals2(ind)-dvals,vals2(ind),Ninterp);
		tmpcpd=interp1(vals,cpprime,tmpval)-cps1(i); %tmpdif=tmpval-vals1(i);
		mini=squeeze(nearest(tmpcpd,p));
		rangedata(i,:)=[vals1(i) median(tmpval(mini)) median(tmpcpd(mini))]; end, end
rangedata=rangedata(~isnan(rangedata(:,1)),:);
% h=fdesign.lowpass('N,F3dB',12,0.15);
% d1=design(h,'butter');
% y=filtfilt(d1.sosMatrix,d1.ScaleValues,medfilt1(diff(rangedata(:,1:2),1,2),3));
%figure; plot(rangedata(:,1),diff(rangedata(:,1:2),1,2),'k.--')
indrange=find(isnear(diff(rangedata(:,1:2),1,2),min(diff(rangedata(:,1:2),1,2))));
valrange=[mean(rangedata(indrange,1)) mean(rangedata(indrange,2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Backward computation %%%
try, %for symmetrical distributions, this forces symmetry in valrange
	if any(ps<0), cps(1)=exp(psprime(end));
		for n=2:length(psprime),
			cps(n)=exp(logsum(str8n(psprime(end-n+1:end)),1,'baseIn','natural')); end
	else cps=cumsum(psprime(end:-1:1)); end
	if any(diff(valsprime,1)~=diff(valsprime(1:2))),
		vals=linspace(valsprime(1),valsprime(end),length(vals)); vals=vals(end:-1:1);
		cpprime=interp1(valsprime(end:-1:1),cps,vals,'pchip'); end
	cpprime=cpprime/cpprime(end);
	vals=-vals;
	
	R1range=[vals(1) vals(find(cpprime>1-p,1))]; vals1=linspace(R1range(1),R1range(2),Ninterp1); cps1=interp1(vals,cpprime,vals1);
	R2range=[vals(find(cpprime<p,1,'last')) vals(end)]; vals2=linspace(R2range(1),R2range(2),Ninterp1); cps2=interp1(vals,cpprime,vals2);
	dcps=str8n(cps2)*ones(1,Ninterp1)-ones(Ninterp1,1)*cps1;
	dvals=diff(vals(1:2));
	
	rangedata=nan(Ninterp1,3);
	for i=1:Ninterp1, ind=find(dcps(:,i)>=p,1);
		if ~isempty(ind),
			tmpval=linspace(vals2(ind)-dvals,vals2(ind),Ninterp);
			tmpcpd=interp1(vals,cpprime,tmpval)-cps1(i);
			mini=nearest(tmpcpd,p);
			rangedata(i,:)=[vals1(i) median(tmpval(mini)) median(tmpcpd(mini))]; end, end
	rangedata=rangedata(~isnan(rangedata(:,1)),:);
	% h=fdesign.lowpass('N,F3dB',12,0.15);
	% d1=design(h,'butter');
	% y=filtfilt(d1.sosMatrix,d1.ScaleValues,medfilt1(diff(rangedata(:,1:2),1,2),3));
	%figure; plot(rangedata(:,1),y,'k.--')
	indrange=find(isnear(diff(rangedata(:,1:2),1,2),min(diff(rangedata(:,1:2),1,2))));
	valrange=(valrange-[mean(rangedata(indrange,2)) mean(rangedata(indrange,1))])/2;
catch
	disp('CIpALERT: forward computation only'); end

%for n=1:2, indrange(n)=round(median(nearest(valsprime,valrange(n)))); end
vals=-vals(end:-1:1);
patborders=interp1(vals,fcps,valrange); end


%BINNER.m
%
% usage: [Bcell Bmat]=binner(Din,Tin,[Tmin Tmax Nbins]); Creates Bcell, Bmat based on [Tmin Tmax Nbins] and populates appropriate bins with Din
%        [Bcell]=binner(Din,Tin,Bcell,Bmat); Populates Bcell with (more, if already populated) Din based in Bmat binnings
%
% Din must be an nx2 dataset with the to-be-mashed data in C1, and the time-values that define binning in C2
% Bcell is an Nx1 cell vector with data in each cell corresponding to the segments of the T axis defined by Bmat
% Bmat is an Nx2 matrix where the lefthand column is the start of each bin in T and the righthand column the end of each bin in T

function [Bcell Bmat]=binner(Din,Tin,Bcell,Bmat)
if ~iscell(Bcell), Bmat=linspace(Bcell(1),Bcell(2),Bcell(3)+1); clear Bcell;
	Bmat=[Bmat(1:end-1)' Bmat(2:end)']; Bcell=cell(size(Bmat,1),1); end

if ~isempty(Din), if length(Din)~=length(Tin), error('data and time vectors must be same length'); end
	if size(Bmat,1)>length(Din),
		for b=1:size(Bmat,1),
			inds=find(and(Bmat(b,1)<Tin,Bmat(b,2)>=Tin));
			Bcell{b}(end+[1:length(inds)])=Din(inds); end
	else
		for d=1:length(Din),
			if d==1, ind=find(and(Bmat(:,1)<=Tin(d),Bmat(:,2)>Tin(d)));
			else ind=find(and(Bmat(:,1)<Tin(d),Bmat(:,2)>=Tin(d))); end
			Bcell{ind}(end+1)=Din(d); end, end, end, end


% vardef.m
%
% Calculate variance of a probability distribution
%
% [V M]=vardef(p,X,plotNY01);
% p: vector of log-probabilities
% X: 
% plotNY01: Activates plotting (OPTIONAL)
%
% TE Hudson
% 10.05.08

function [V M]=vardef(Lp,X,Fnum,spacing)
if nargin<4, spacing=1; end
if nargin<3, Fnum=[]; end
X=X(:); if length(X)==2, X=linspace(X(1),X(2),size(Lp,1))'; end
if and(size(Lp,1)==1,size(Lp,2)>1), Lp=Lp'; end

for c=1:size(Lp,2), if or(isnear(sum(Lp(:,c)),1),isnear(max(Lp(:,c)),1)), Lp(:,c)=log(Lp(:,c)/nansum(Lp(:,c))); %normalize and log
	else inds=find(isnear(Lp(:,c),max(Lp(:,c))));
		if isnan(max(Lp(:,c))), error('input likelihoods are nan'); end
		if or(any(inds(1)==[1 length(X)]),any(inds(end)==[1 length(X)])), Lp(:,c)=-Lp(:,c); end %reverse if input logP is upside-down
		Lp(:,c)=Lp(:,c)-logsum(Lp(:,c)); end, end %normalize the log

%M=meandef(Lp,X,Fnum,spacing)
M=logpeakval([Lp X]);

if ~all(isnear(diff(X),diff(X(1:2)))), %if xs are not linearly spaced, re-space xs and interpolate ps
	Xls=linspace(X(1),X(end),length(X))';
	Lpls=interp1(X,Lp,Xls); Lpls=Lpls-logsum(Lpls);
else Xls=X; Lpls=Lp; end
LX2=log((Xls*ones(1,size(Lpls,2))-ones(size(Lpls,1),1)*M).^2);
V=exp(logsum(LX2*ones(1,size(Lpls,2))+Lpls,1,spacing,'baseIn','natural'));

if ~isempty(Fnum), if Fnum==0, figure; else figure(Fnum); end
	hold on;
	for n=1:size(Lp,2),
		plot(X,Lp(:,n),'k-');
		plot(M(n)*[1 1],[min(Lp(:,n)) min(Lp(:,n))/4],'r-')
		plot(M(n)+sqrt(V(n))*[-1 1],min(Lp(:,n))/4*[1 1],'r-');
	end, r=axis; axis([r(1) r(2) min(min(Lp)) 0]); end, end


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
	else plot(inds,meanvec(inds),'ko--','MarkerFaceColor',COLface,'MarkerEdgeColor',COLedge,'MarkerSize',7); end, end, end


function tf=isnear(a,b,tol)
%ISNEAR True Where Nearly Equal.
% ISNEAR(A,B) returns a logical array the same size as A and B that is True
% where A and B are almost equal to each other and False where they are not.
% A and B must be the same size or one can be a scalar.
% ISNEAR(A,B,TOL) uses the tolerance TOL to determine nearness. In this
% case, TOL can be a scalar or an array the same size as A and B.
%
% When TOL is not provided, TOL = SQRT(eps).
%
% Use this function instead of A==B when A and B contain noninteger values.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% Mastering MATLAB 7
% 2005-03-09

%--------------------------------------------------------------------------
if nargin==2
   tol=sqrt(eps);
end
if ~isnumeric(a) || isempty(a) || ~isnumeric(b) || isempty(b) ||...
   ~isnumeric(tol) || isempty(tol)
   error('Inputs Must be Numeric.')
end
if any(size(a)~=size(b)) && numel(a)>1 && numel(b)>1
   error('A and B Must be the Same Size or Either can be a Scalar.')
end
tf=abs((a-b))<=abs(tol); end


%function straightens a matrix (useful in nested arguments)
function [vecout]=str8n(matin)
vecout=matin(:); end





