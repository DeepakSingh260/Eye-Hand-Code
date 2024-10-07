function Scal=calFinger
global snd Sexp Tx wndw

load C:\Users\VMILlab\Desktop\polhemusCalData.mat Scal P

clc
DrawFormattedText(wndw,'Finger Calibration','center',Sexp.resy/3,Sexp.col.BLU); Screen('Flip', wndw); 
fprintf('\n'); fprintf('\n'); fprintf('\n'); fprintf('\n');
cprintf('*RED','Instruct subject to place finger on the button','(press any key)')
KbStrokeWait; count=1; nogood=1;
while KbCheck; end %until key is pressed
while nogood, clc
	fprintf('\n'); fprintf('\n'); 
	cprintf('*BLUE',['Attempt #' num2str(count) '\n']); fprintf('\n');
	cprintf('red','collecting...\n');
    play(snd.CLICK); pause(.3); play(snd.CLICK); pause(.3); play(snd.CLICK); pause(1);
	
    sampmat=pollPHstream(Scal.Clist);
	sampmat(end+1:Scal.Nsamp,:)=nan;
	for n=2:Scal.Nsamp,
		sampmat(n,:)=pollPHstream(Scal.Clist); end
    Scal.mS=mean(sampmat,1)';
	Scal.AE=(Scal.T*Scal.mS-Scal.F')';
    Scal.AAE=mean(abs((Scal.T*sampmat')'-ones(Scal.Nsamp,1)*Scal.F));
    Scal.SD=std((Scal.T*sampmat')',0,1);
    nogood=any(Scal.SD>1);
    
	fprintf('\n'); fprintf('\n'); 
	cprintf('*BLUE',['Attempt #' num2str(count)]); 
	if nogood, count=count+1; snd.TIMEOUTPTHH(1);
		cprintf('*RED',' -- FAILURE!\n'); fprintf('\n');
	else
		%Scal.scale=[Sexp.xy0./mean(Scal.Dxy,1) 1];
        Scal.scale=[Sexp.xy0./(Sexp.MMxy/2) 1];
		Tx=@(newsamp) Scal.scale.*([Scal.T Scal.F'-Scal.T*Scal.mS-Scal.C']*[newsamp 1]')';
		cprintf('*RED',' -- SUCCESS!!\n'); fprintf('\n'); end
	cprintf('blue',['Avg Err=[' num2str(round(Scal.AE(1)*100)/100) ', ' num2str(round(Scal.AE(2)*100)/100) ', ' num2str(round(Scal.AE(3)*100)/100) '] \n']);
	cprintf('blue',['Avg AbsErr=[' num2str(round(Scal.AAE(1)*10)/10) ', ' num2str(round(Scal.AAE(2)*10)/10) ', ' num2str(round(Scal.AAE(3)*10)/10) '] \n']);
	cprintf('blue',['SD=[' num2str(round(Scal.SD(1)*100)/100) ', ' num2str(round(Scal.SD(2)*100)/100) ', ' num2str(round(Scal.SD(3)*100)/100) '] \n']);
	
	if ~nogood,
		fprintf('\n'); fprintf('\n');
		cprintf('*blue','As a final check, the following should be near the origin: \n')
		cprintf('*black',['    P=[' num2str(Tx(P(5,Scal.Clist))) '] \n']); end
	pause(4.1); end
Sexp.Clist=Scal.Clist;
save C:\Users\VMILlab\Desktop\polhemusCalData.mat sampmat Scal -APPEND;