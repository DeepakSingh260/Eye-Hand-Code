
function tstVG(Spre,el) %#ok<*NASGU,*AGROW>
global Sexp ELxy PHxy Ftmp wndw snd
Sexp=Spre; Sexp.quitExp=0;
Screen('Flip', wndw);

Nprob=[0 .5; .5 .79; .79 .93; .93 1];
snd.TINK=audioplayer(15*sin(1.2*pi*0.05*[0:50]),20000);

%%%%%%%%%%%%%%%%%%%%%%%
%%% SAC Sexp params %%%
Sexp.ReachDist=101.5; %mm (4in)
Sexp.Ntarg=5;
Sexp.Lc=1;
if strcmp(Sexp.PhaseNow,'S1'),
    Sexp.Nsac=25;
    Sexp.OArot=0;
elseif strcmp(Sexp.PhaseNow,'S2'),
    Sexp.Nsac=25;
    Sexp.OArot=pi/4;
elseif strcmp(Sexp.PhaseNow,'S3'),
    Sexp.Nsac=25;
    Sexp.OArot=pi/7; end
Sexp.BTlist=zeros(Sexp.Nsac,1);
Sexp.EXPlist=Shuffle(str8n([1:Sexp.Ntarg]'*ones(1,Sexp.Nsac/Sexp.Ntarg)));
Sexp.Tlist=rot(rot([Sexp.ReachDist 0],2*pi/Sexp.Ntarg*[0:Sexp.Ntarg-1]),Sexp.OArot); %mm
EXP=zeros(Sexp.Nsac,22); %startXY(1:2) targXY(3:4) startEL1(5:6) startEL2(7:8) startPH(9:10) endPH(11:12) endEL1(13:14) endEL2(15:16) tstartPH(17) tendPH(18) tstartEL1(19) tendEL1(20) tstartEL2(21) tendEL2(22)
EXP(:,3:4)=Sexp.Tlist(Sexp.EXPlist,:); %mm

Sexp.Tfixup=cell(Sexp.Nsac,1);
Sexp.Tgo=nan;
Sexp.ELgo=nan;
Sexp.sampintervalsec=7;
Sexp.keypress=cell(Sexp.Nsac,1);
Sexp.isacstart=nan(Sexp.Nsac,2); Sexp.ifinstart=nan(Sexp.Nsac,1);
Sexp.isacend=nan(Sexp.Nsac,2); Sexp.ifinend=nan(Sexp.Nsac,1);

%Fixation params
sz=round(.9*el.calibrationtargetsize/100*Sexp.resx);
inset=round(.85*el.calibrationtargetwidth/100*Sexp.resx);

%Fingertip Indicator
Frect=CenterRectOnPoint([0 0 21 21],0,0);
Fposs=[0 10 10 -10 -10; 0 10 -10 10 -10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Eye Position Spotlight %%%
si=105; %size
tw=4*si+1; th=4*si+1; %Size of support in pixels, derived from si:
contrast=2.5; %Contrast
aspectratio=1; %Aspect ratio
sc=45; %Spatial constant of the exponential "hull"
blobtex=CreateProceduralGaussBlob(wndw,tw,th);
texrect=Screen('Rect',blobtex);
myblobpars=[contrast,sc,aspectratio,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saccade Start Screen %%%
if strcmp(Sexp.PhaseNow,'S3'),
    DrawFormattedText(wndw,'Touch targets while fixating screen center','center',Sexp.resy/3,Sexp.col.RED); Screen('Flip', wndw);
else DrawFormattedText(wndw,'Look at & touch the dot cloud','center',Sexp.resy/3,Sexp.col.BLU); Screen('Flip', wndw); end
clc; cprintf('*blue',['press any button to begin saccade phase #' Sexp.PhaseNow '\n']); KbStrokeWait;

%MACkeycodes: shift: 225, ctrl: 224, tab: 43, capslock: 57
%PC keycodes: shift: 16,  ctrl: 17,  tab: 9,  capslock: 20
keys.shift=16; keys.ctrl=17; keys.tab=9; keys.capslock=20; KbCheck; %#ok<STRNU> %clear key buffer
Sexp.keys.backspaceCode=8;

%%%%%%%%%%%%%%%%
%%% MAIN EXP %%%
%%%%%%%%%%%%%%%%
t=0; BTnow=0; ELdat=cell(0); PHdat=cell(0); Sexp.tSAC=GetSecs;
while t<Sexp.Nsac, t=t+(~BTnow); BTnow=0; %reset bad-trial flag for next trial

    [keydown, ~, keyCode]=KbCheck; %#ok<ASGLU>
	Sexp.Fnow=EXP(t,1:2); %mm
	if isfield(Sexp,'Tnow'), Sexp=rmfield(Sexp,'Tnow'); end
	Fpix=Sexp.mm2pix(Sexp.Fnow);
	startargetrectO=CenterRectOnPoint([0 0 sz sz],Fpix(1),Fpix(2));
	startargetrectI=CenterRectOnPoint([0 0 .3*[sz sz]],Fpix(1),Fpix(2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Gaussian dot-cloud %%%
	sigma=31*[1 0;0 1]; Ndots=65; dtSmooth=30; dtG=.06; countermax=round(300/dtG)+1;
	Dxy=zeros(Ndots,2,countermax); Tpix=Sexp.mm2pix(EXP(t,3:4)); Lthresh=dist(Tpix,Fpix,2)/3; Vthresh=Sexp.deg2pix([75 0]);
	for n=1:round(dtSmooth/dtG)+1,
		Dxy(:,:,n)=mvnrnd(Sexp.mm2pix(EXP(t,3:4)),sigma,Ndots); end
	SZdxy3=size(Dxy,3);
	tlistDxy=linspace(0,Sexp.sampintervalsec,round(Sexp.sampintervalsec/dtG)); tlistDxy=[tlistDxy tlistDxy(end)];
	Sexp.Scol=[.7*Sexp.col.WHT; Sexp.col.BLK; Sexp.col.GRY; Sexp.col.BLK]; Sexp.Sr=[13; 13-.7; 7.5; 2.5];
	counter=1;
	    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FIXATION bet 250-750 ms %%%%
    keyDown=0;
    while ~keyDown,
        clc; cprintf('cyan','fixation loop...\n');
        cprintf('blue','press L-arrow to drift correct \n');
        cprintf('blue','press R-arrow to recalibrate \n');
        cprintf('blue',['trial # ' num2str(t)  ' \n']);
        
        sampinitPH(Sexp.sampintervalsec); sampinitEL(Sexp.sampintervalsec); %pre-allocate memory for sampintervalsec of data
        Sexp.Sdelay(t)= 1; %const delay
        Sexp.Tfixup{t}=GetSecs; %first entry
        pollPH; pollEL;
        while and(~keyDown,GetSecs-Sexp.Tfixup{t}(end)<Sexp.Sdelay(t)), %Wait for delay (keyDown event breaks loop)
            Screen('DrawTextures',wndw,blobtex,[],CenterRectOnPoint(texrect,ELxy(1),ELxy(2))',0,[],[],[],[],kPsychDontDoRotation,myblobpars');
            Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
            Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
            Screen('FillOval',wndw,Sexp.col.RED,Frect+[PHxy(1) PHxy(2) PHxy(1) PHxy(2)]); %show Ftip
            pollPH; Screen('Flip',wndw);
            [keyDown,~,kbNow]=KbCheck;
            
            %%% PRE-CONDITIONING so that Ftip in particular is stationary when target comes up
            othercrit=1;
            if Sexp.jF>Sexp.NthreshPH,
                if std(dist(Sexp.pix2mm(Ftmp(Sexp.jF+[-Sexp.NthreshPH:0],1:2)),EXP(t,1:2),2))>Sexp.FprethreshSD,
                    othercrit=0; end, end
            if or(~othercrit,or(dist(ELxy(1:2),Sexp.xy0,2)>55,dist(PHxy(1:2),Sexp.xy0,2)>30)),
                Sexp.jE=0; Sexp.jF=0;
                Sexp.Tfixup{t}(end+1)=GetSecs; end %add new start time because finger is not reset properly
            pollEL; pollPH; end
        
        
        Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
        Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
      
        
        
        Screen('Flip',wndw);
        
        %%%%%%%%%%%%%%%%%%%%%
        %%% Start New Saccade(SPACE) or do drift correct(R-arrow key) or recalibrate(L-arrow key)
        if keyDown, %do drift correct or recalibration
            if find(kbNow)==Sexp.keys.lrCode(1), %drift correction (left-arrow)
                clear KbCheck
                Eyelink('StopRecording');
                clc; fprintf('\n'); cprintf('red','press SPACEBAR to perform Drift Correction\n');
                Eyelink('DriftCorrStart',Sexp.x0,Sexp.y0,1,1,0);
                EyelinkClearCalDisplay(el);
                [~, DCmessage]=Eyelink('CalMessage'); Sexp.DCdat{t,1}=DCmessage;
                Eyelink('ApplyDriftCorr');
                Eyelink('StartRecording');
                keyDown=0; %in either case, re-initiate the variable fixation period
    
            elseif find(kbNow)==Sexp.keys.lrCode(2), %re-calibrate (right-arrow)
                doELcal(el);  keyDown=0; %in either case, re-initiate the variable fixation period
            elseif find(kbNow)==Sexp.keys.backspaceCode, keyDown=1; Sexp.quitExp=1; t=Sexp.Nsac; end %quit exp
        else keyDown=1; end, end %begin saccade sequence
    clc; cprintf('blue',['trial # ' num2str(t)  ' \n']);
    
    if ~Sexp.quitExp,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Give countdown to go %%%
    onnow=0; Sexp.t0(t)=GetSecs; Sexp.jE=0; Sexp.jF=0; pollEL; pollPH; %start of recording session
    for n=1:3,
        while GetSecs-Sexp.t0(t)<n/2,
            if and(GetSecs-Sexp.t0(t)>.65,GetSecs-Sexp.t0(t)<1.15), onnow=1;
                if GetSecs-Sexp.t0(t)>tlistDxy(counter), counter=(counter~=SZdxy3)*counter+1;
                    Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
                    Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
                    Screen('DrawDots',wndw,Dxy(:,:,counter)',3,Sexp.col.BLU,[],2); %target cloud
                    Screen('Flip',wndw); end
            elseif onnow, onnow=0;
                Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
                Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
                Screen('Flip',wndw); end
            pollPH; end, end
    if strcmp(Sexp.PhaseNow,'S1'), playblocking(snd.TINK); end
    Screen('Flip',wndw);
    
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%% LOOK & REACH %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Sexp.Tnow=EXP(t,3:4); 
    goodindex=[0 0 0];
    Sexp.ithreshF=[0 0]; %reset start/end indices for next trial
    startz=PHxy(3); Sexp.Tgo(t)=GetSecs; Sexp.ELgo(t)=GetSecs;

    
    if strcmp(Sexp.PhaseNow,'S2'),    
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% Perform Saccades %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Sexp.Tsac{t}=zeros(0,2); pnow=[1 1]'*[Fpix GetSecs];
        % Nsac=0; randnow=rand; Sexp.Scount(t)=2*(find(and(Nprob(:,1)<randnow,Nprob(:,2)>randnow))-1);
        Sexp.Sdt(t)=1.5+rand; counter=1; tcounter=GetSecs;
        while GetSecs-Sexp.Tgo(t)<Sexp.Sdt(t), pollPH; 
            
            %%%%%%%%%%%%%%%%%%%%%%%
            %%% Maintain Target %%%
            if GetSecs-tcounter>tlistDxy(counter), counter=counter+1;
                if counter==countermax, counter=1; tcounter=GetSecs; end
                Screen('DrawDots',wndw,Dxy(:,:,counter)',3,Sexp.col.BLU,[],2); %target cloud
                Screen('Flip',wndw); end; end, end
    pollEL
 
    
    
    if ~strcmp(Sexp.PhaseNow,'S1'),
        play(snd.TINK); end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% FINGER FIRST %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Sexp.Tbeep(t)=GetSecs; counter=1; tcounter=Sexp.Tbeep(t);
    while GetSecs-Sexp.Tbeep<8, %or(Sexp.ithreshF(2)==0,Sexp.jF<Sexp.ithreshF(2)),
        %if GetSecs-Sexp.Tbeep(t)>Sexp.ELtimeout+(strcmp(Sexp.Stype,'C')), BTnow=1; %3-second rule
        %    clc; cprintf('red','3-second rule -- trial incomplete \n'); break; end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GET DATA AND UPDATE SCREEN IMAGE %%%
        pollPH; pollEL;
        Screen('DrawTextures',wndw,blobtex,[],CenterRectOnPoint(texrect,ELxy(1),ELxy(2))',0,[],[],[],[],kPsychDontDoRotation,myblobpars');
        Screen('FillOval',wndw,Sexp.col.RED,Frect+[PHxy(1) PHxy(2) PHxy(1) PHxy(2)]); %show Ftip
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% WAITING for the finger mvmt to begin & end %%%
        %stdcheck(0); %set Sexp.ithreshF as Ftip mvmt starts and ends
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Maintain Target %%%
        if counter>=117, counter=1; tcounter=GetSecs; end
        if GetSecs-tcounter>tlistDxy(counter), counter=counter+1;
            if counter==countermax-1, counter=1; tcounter=GetSecs; end
            Screen('DrawDots',wndw,Dxy(:,:,counter)',3,Sexp.col.BLU,[],2); %target cloud
            Screen('Flip',wndw); end; end, end, end
%datSAVE(ELdat,PHdat,EXP);
