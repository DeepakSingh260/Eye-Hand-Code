%lisstest.m
% inputs are ELdat(t,:) and PHdat{t}
% output adds a row to the LAV and LPV fields of the Sexp struct

function lisstest(ELdat,PHdat,t)
global Sexp

Fpos=Sexp.Fpos{t}(1:end-1,:);
tF=Fpos(:,3);

%Pcell=cell(1,2); Ecell=cell(1,2);
Pdat=PHdat(:,[1 2]);
Pdat(:,1:2)=rot(Pdat(:,1:2),Sexp.OArot(t));
Pdat(:,1:2)=[Pdat(:,1)*Sexp.SGN(t,1) Pdat(:,2)*Sexp.SGN(t,2)];
tP=PHdat(:,4)-PHdat(1,4);
%dtP=tP(end)-tP(1);

Edat=zeros(size(ELdat{Sexp.EYEnow(1)},1),2);
for Enow=Sexp.EYEnow, Edat=Edat+ELdat{Enow}(:,1:2)/Sexp.Neye; end
Edat(:,1:2)=rot(Edat(:,1:2),Sexp.OArot(t));
Edat(:,1:2)=[Edat(:,1)*Sexp.SGN(t,1) Edat(:,2)*Sexp.SGN(t,2)];
tE=ELdat{Sexp.EYEnow(1)}(:,4)-ELdat{Sexp.EYEnow(1)}(1,4);
%dtE=tE(end)-tE(1);
%Vcell=cell(3,1); 
%Vcell{3}=linspace(-pi/2,pi/2,61);

%for XY=1:2,
   
    
    
    %tmpy=interp1(PHdat{2*n}(:,4)-PHdat{2*n}(1,4),PHdat{2*n}(:,2),Sexp.Fpos{n}(1:end-1,3)); dy=tmpy-Sexp.Fpos{n}(1:end-1,2);
    %tmpx=interp1(PHdat{2*n}(:,4)-PHdat{2*n}(1,4),PHdat{2*n}(:,1),Sexp.Fpos{n}(1:end-1,3)); dx=tmpx-Sexp.Fpos{n}(1:end-1,1);
   
    %hist(dist([dx(:) dy(:)],2),length(dy)/20);
    %Sexp.lissvar(t)=median(dist([dx(:) dy(:)],2))/20;
    %r=axis; plot(Sexp.lisstest(t)*[1 1],r(3:4),'c--','LineWidth',2); end
    
%     
%     
%     
%     
%     
%     
%     
%     Vcell{1}=Sexp.fXY(t,XY)+linspace(-.1,.1,31);
%     Vcell{2}=Sexp.Aprime(t,XY)*linspace(0,4,61);
%     Sexp.FAPdat{t,XY}=Vcell;
%     
%     Pcell{2}=[Pdat(:,XY) tP/dtP];
%     Sexp.FAPdatP{t,XY}=Pcell;
%     try [FAPhat]=Lfap(Pcell,Vcell,2,'Fnum',[],'Fprime',Sexp.fXY(t,XY),'Aprime',Sexp.Aprime(t,XY),'Pprime',0);
%         %inds=[findnearestN(Lcell{3}(:,2),-pi/2,1) findnearestN(Lcell{3}(:,2),pi/2,1)];
%         Sexp.LAV(t,XY)=(FAPhat.amp.gain(1)-1)*diff(FAPhat.amp.gain(end-1:end));
%         Sexp.LPV(t,XY)=FAPhat.phase.rad(1)*diff(FAPhat.phase.rad(end-1:end)); %logpeakval(Lcell{3}(inds(1):inds(2),:))*diff(FAPhat.phase.rad(end-1:end));
%         Sexp.FAPhatP{t,XY}=FAPhat;
%     catch, Sexp.LAV(t,XY)=10; Sexp.LPV(t,XY)=10; Sexp.FAPhatE{t,XY}=nan; end
%     
%     
%     %ERROR MEASURE 1
%     tmpx=interp1(PHdat{2*n}(:,4)-PHdat{2*n}(1,4),PHdat{2*n}(:,1),Sexp.Fpos{n}(1:end-1,3)); dx=tmpx-Sexp.Fpos{n}(1:end-1,1);
%     tmpy=interp1(PHdat{2*n}(:,4)-PHdat{2*n}(1,4),PHdat{2*n}(:,2),Sexp.Fpos{n}(1:end-1,3)); dy=tmpy-Sexp.Fpos{n}(1:end-1:2);
%     Sexp.Dxy(t,XY)=mean(dist([dx dy],Sexp.Fpos{n}(1:end-1,1:2),2));
% 
%     %ERROR MEASURE 2
%     Ecell{2}=[Edat(:,XY) tE/dtE];
%     Sexp.FAPdatE{t,XY}=Ecell;    
%     try [FAPhat]=Lfap(Ecell,Vcell,2,'Fnum',[],'Fprime',Sexp.fXY(t,XY),'Aprime',Sexp.Aprime(t,XY),'Pprime',0);
%         %inds=[findnearestN(Lcell{3}(:,2),-pi/2,1) findnearestN(Lcell{3}(:,2),pi/2,1)];
%         Sexp.LAV(t,XY+2)=(FAPhat.amp.gain(1)-1)/diff(FAPhat.amp.gain(end-1:end));
%         Sexp.LPV(t,XY+2)=FAPhat.phase.rad(1)/diff(FAPhat.phase.rad(end-1:end)); %logpeakval(Lcell{3}(inds(1):inds(2),:))*diff(FAPhat.phase.rad(end-1:end));
%         Sexp.FAPhatE{t,XY}=FAPhat;
%     catch, Sexp.LAV(t,XY+2)=10; Sexp.LPV(t,XY+2)=10; Sexp.FAPhatE{t,XY}=nan; end
% Sexp.lissvar(t)=10*dist([Sexp.LAV(t,:) Sexp.LPV(t,:)/pi],2); 

 %end
 
%EYE ERROR
tmpxE=interp1(tE,medfilt1(Edat(:,1),15),tF); dxE=tmpxE-Fpos(:,1);
tmpyE=interp1(tE,medfilt1(Edat(:,2),15),tF); dyE=tmpyE-Fpos(:,2);
Sexp.lissvar(t,1)=nanstd(dist([dxE(:) dyE(:)],2)/25);

%FINGER ERROR
tmpxF=interp1(tP,medfilt1(Pdat(:,1),7),tF); dxF=tmpxF-Fpos(:,1);
tmpyF=interp1(tP,medfilt1(Pdat(:,2),7),tF); dyF=tmpyF-Fpos(:,2);
Sexp.lissvar(t,2)=nanstd(dist([dxF dyF],2)/25);

save tmpdatLV


% %%%%%%%%%%%%%%%%%%%%%%
% %%% testing script %%%
% %  first load Sexp  %
% ploton=1;
% if ploton, Fnum{1}=17; Fnum{2}=19; else Fnum{1}=[]; Fnum{2}=[]; end
% for t=1:length(Sexp.lissvar),
%     for XY=1:2,
%         [FAPhat Lcell]=Lfap(Sexp.FAPdatP{t,XY},Sexp.FAPdat{t,XY},2,'Fnum',Fnum{1},'Fprime',Sexp.fXY(t,XY),'Aprime',Sexp.Aprime(t,XY),'Pprime',0);
%         inds=[findnearestN(Lcell{3}(:,2),-pi/2,1) findnearestN(Lcell{3}(:,2),pi/2,1)];
%         Sexp.LAV(t,XY)=(FAPhat.amp.gain(1)-1)*diff(FAPhat.amp.gain(end-1:end));
%         Sexp.LPV(t,XY)=logpeakval(Lcell{3}(inds(1):inds(2),:))*diff(FAPhat.phase.rad(end-1:end));
%         Sexp.FAPhatP{t,XY,1}=FAPhat;
%         
%         [FAPhat Lcell]=Lfap(Sexp.FAPdatE{t,XY},Sexp.FAPdat{t,XY},2,'Fnum',Fnum{2},'Fprime',Sexp.fXY(t,XY),'Aprime',Sexp.Aprime(t,XY),'Pprime',0);
%         inds=[findnearestN(Lcell{3}(:,2),-pi/2,1) findnearestN(Lcell{3}(:,2),pi/2,1)];
%         Sexp.LAV(t,2+XY)=(FAPhat.amp.gain(1)-1)/diff(FAPhat.amp.gain(end-1:end));
%         Sexp.LPV(t,2+XY)=logpeakval(Lcell{3}(inds(1):inds(2),:))/diff(FAPhat.phase.rad(end-1:end));
%         Sexp.FAPhatE{t,XY,2}=FAPhat; end
%     Sexp.lissvar(t)=10*sqrt(sum([Sexp.LAV(t,:) Sexp.LPV(t,:)/pi].^2)); end
% D=nan(2); d=nan(2); S=nan(2);
% Dxy=[];
% for n=1:2,
% tmpx=interp1(PHdat{2*n}(:,4)-PHdat{2*n}(1,4),PHdat{2*n}(:,1),Sexp.Fpos{n}(1:end-1,3)); dx=tmpx-Sexp.Fpos{n}(1:end-1,1);
% tmpy=interp1(PHdat{2*n}(:,4)-PHdat{2*n}(1,4),PHdat{2*n}(:,2),Sexp.Fpos{n}(1:end-1,3)); dy=tmpy-Sexp.Fpos{n}(1:end-1,2); 
% Dxy(n)=mean(dist([dx dy],Sexp.Fpos{n}(1:end-1,1:2),2));
% 
% d(n,:)=[mean(dx) mean(dy)];
% D(n,:)=[mean(abs(dx)) mean(abs(dy))];
% S(n,:)=sqrt([var(dx) var(dy)]);
% end



