function [Scal]=calTable(forceCal)
if nargin==0, forceCal=0; end
global snd Sexp wndw

% init Sexp fields
makeSexp;

% Define sounds
snd.GOBEEP=audioplayer(sin(2*pi*.1*[0:1000]),15000);
snd.GOODBEEP=audioplayer([sin(2*pi*.07*[0:3:2000]) sin(2*pi*.12*[0:3:4000])],8000);
snd.RASPBERRY=audioplayer(sin(2*pi*0.004*[0:2:8000])*10,8000);
snd.TIMEOUTPTHH=audioplayer([sin(2*pi*0.03*[0:3:500]) sin(2*pi*0.028*[0:3:500]) sin(2*pi*0.026*[0:3:500]) sin(2*pi*0.024*[0:3:500]) sin(2*pi*0.022*[0:3:500]) sin(2*pi*0.02*[0:3:500]) sin(2*pi*0.018*[0:3:500]) sin(2*pi*0.016*[0:3:500]) sin(2*pi*0.014*[0:3:500]) sin(2*pi*0.012*[0:3:500]) sin(2*pi*0.01*[0:3:500]) sin(2*pi*0.008*[0:3:500]) sin(2*pi*0.006*[0:3:500]) sin(2*pi*0.004*[0:3:500])],8000);
snd.CLICK=audioplayer(sin(1.4*pi*0.05*[0:200]),45000); play(snd.CLICK); %initialize PsychPortAudio with first snd command

%
doCal=1;
 % if ~forceCal,
 % 	if exist('C:\Users\VMILlab\Desktop\polhemusCalData.mat','file'), load C:\Users\VMILlab\Desktop\polhemusCalData.mat Scal lastcaldate %lastssid
 % 		if datenum(lastcaldate)==datenum(date), doCal=0; end, end, end

if doCal, Scal.calgood=0; lastcaldate=date;
	while ~Scal.calgood,
	sampmove=cell(8,1); lastcaldate=date; %RUN THE CAL PROCEDURE
	Scal.Cloc=[[206.375*[-1 0 1 -1 0 1 -1 0 1] -293.6875]' [119.0625*[1 1 1 0 0 0 -1 -1 -1] -147.6375]'];
	Tr=[15;10;5]; Tcol1=[Sexp.col.GRY/2; Sexp.col.BLU; Sexp.col.GRY/2]; Tcol2=[Sexp.col.GRY/2; Sexp.col.RED; Sexp.col.GRY/2];
	
    P=nan(10,Sexp.nPHsamp);
	for p=1:10, %drawTarg(Scal.Cloc(1:9,:),5,Sexp.col.BLU,0); %target grid
		drawTarg(Scal.Cloc(p,:),Tr,Tcol2,1);   %current target
		play(snd.GOBEEP); %movement go-signal for next location
					
		
		%FIRST TARGET
		if p==1, sampmat=nan(Sexp.nPHsamp);
            Screen('Flip', wndw); 
            DrawFormattedText(wndw,'Table Calibration','center','center',Sexp.col.BLU); 
			drawTarg(Scal.Cloc(p,:),Tr,Tcol2,1);
            cprintf('RED*','press key when first position is covered and finger is still\n');
			while ~KbCheck, %until key is pressed
				sampmat=[pollPHstream(1:Sexp.nPHsamp); sampmat];
                
            sampmat(1,2:6)%-sampmat(2,5)
            end
           
       
     
            
            tnowx=GetSecs; sampmatTF=sampmat(:,Sexp.nPHsamp)>GetSecs-1.1;
            save SMAT.mat sampmat Sexp tnowx sampmatTF
            
			sampmat=sampmat(sampmat(:,Sexp.nPHsamp)>GetSecs-1.1,:); %all data for last 1.1s
			Scal.vars=1.25*var(sampmat(:,1:end-1)); Scal.Nsamp=size(sampmat,1); Scal.Csamp=size(sampmat,2)-1;
		    varmult=(Scal.vars>0); clc
			
		%ALL OTHER TARGETS	
		else
			%wait for movement
			sampmat=nan(0,Scal.Csamp+1);
			while size(sampmat,1)<Scal.Nsamp, %fill matrix
				sampnow=pollPHstream(1:Sexp.nPHsamp);
				sampmat=[sampnow; sampmat]; end
			%check for var thresh to pass (movement ongoing)
			while and(~KbCheck,all(varmult.*var(sampmat(end-Scal.Nsamp+1:end,1:end-1),1,'omitnan')<=2*Scal.vars)), %if any var increases beyond 1.25*Scal.vars
				sampnow=pollPHstream(1:Sexp.nPHsamp);
                disp([max(var(sampmat(end-Scal.Nsamp+1:end,1:end-1),1,'omitnan')) 2*Scal.vars])
				sampmat=[sampnow; sampmat(1:Scal.Nsamp,:)]; end
			sampmove{p}=sampmat; sampmat=nan(Scal.Nsamp,Scal.Csamp+1);
			drawTarg(Scal.Cloc(p,:),Tr,Tcol1,1); play(snd.CLICK);
			
			%hold sensor still for location measurement
			while and(~KbCheck,or(isnan(sampmat(end,1)),any(2*varmult.*var(sampmat(:,1:end-1),1,'omitnan')>Scal.vars))),
				sampmat(2:end,:)=sampmat(1:end-1,:);
				sampmat(1,:)=pollPHstream(1:Sexp.nPHsamp); end, end
		P(p,:)=nanmedian(sampmat,1); Pcell{p}=sampmat; disp('added to pcell');
		play(snd.GOODBEEP); disp('played goodbeep'); end
    
    disp('passed daq-loop'); 
	[Scal]=initStream(Pcell,sampmove,Scal);    
    save C:\Users\VMILlab\Desktop\polhemusCalTemp.mat Scal lastcaldate P Pcell;
    
    
    foo
    end %THIS ONE ADDS THE Z-DIM to the XY, defining Scal.Clist as the triple of xyz-columns in the raw data stream
delete C:\Users\VMILlab\Desktop\polhemusCalTemp.mat
save C:\Users\VMILlab\Desktop\polhemusCalData.mat Scal lastcaldate P Pcell;

end
%lastssid=Sexp.Sinit;
%save zubs\ztmp\polhemusCalData.mat lastssid -APPEND;
