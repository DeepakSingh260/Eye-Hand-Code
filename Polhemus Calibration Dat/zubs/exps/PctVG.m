
function Spre=PctVG(Spre) %#ok<*NASGU,*AGROW>
global Sexp PHxy Ftmp wndw snd
Sexp=Spre; Sexp.quitExp=0;
Screen('Flip', wndw);

Nprob=[0 .5; .5 .79; .79 .93; .93 1];
snd.TINK=audioplayer(15*sin(1.2*pi*0.05*[0:50]),20000);

%%%%%%%%%%%%%%%%%%%%%%%
%%% SAC Sexp params %%%
Sexp.ReachDist=101.5; %mm (4in)
Sexp.Ntarg=10;
Sexp.Lc=1;
Sexp.Nsac=2*Sexp.Ntarg; %Nsac=5*Sexp.Ntarg; %must be multiple of Ntarg
Sexp.OArot=-pi/2;
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

%Fingertip Indicator
Frect=CenterRectOnPoint([0 0 21 21],0,0);
Fposs=[0 10 10 -10 -10; 0 10 -10 10 -10];
%Fixation params
sz=2*3*Fposs(1,2);

%%%%%%%%%%%%%%%%%%%%
%%% Start Screen %%%
DrawFormattedText(wndw,'Touch the dot cloud','center',Sexp.resy/3,Sexp.col.BLU); Screen('Flip', wndw);
clc; cprintf('*blue',['press any button to begin practice']); KbStrokeWait;

%MACkeycodes: shift: 225, ctrl: 224, tab: 43, capslock: 57
%PC keycodes: shift: 16,  ctrl: 17,  tab: 9,  capslock: 20
keys.shift=16; keys.ctrl=17; keys.tab=9; keys.capslock=20; KbCheck; %#ok<STRNU> %clear key buffer
Sexp.keys.backspaceCode=8;
blbCOL=2*Sexp.col.RED/3;
blbCOL(2,:)=round(1.1*mean(blbCOL));

%%%%%%%%%%%%%%%%
%%% MAIN EXP %%%
%%%%%%%%%%%%%%%%
t=0; BTnow=0; ELdat=cell(0); PHdat=cell(0); Sexp.tSAC=GetSecs; Sexp.rTpix=[]; rTpix=75; Sexp.NrT=7; Sexp.Fend=rot([0 rTpix/sqrt(2)],2*pi*rand(Sexp.NrT,1));
while t<Sexp.Nsac, 
    t=t+(~BTnow); BTnow=0; %reset bad-trial flag for next trial
    
    keydown=KbCheck;
    Sexp.Fnow=EXP(t,1:2); %mm
    if isfield(Sexp,'Tnow'), Sexp=rmfield(Sexp,'Tnow'); end
    Fpix=Sexp.mm2pix(Sexp.Fnow);
    startargetrectO=CenterRectOnPoint([0 0 sz sz],Fpix(1),Fpix(2));
    startargetrectI=CenterRectOnPoint([0 0 .8*[sz sz]],Fpix(1),Fpix(2));
    
    %%%%%%%%%%%%%%%%%%%%%
    %%% Feedback Blur %%%
    Sexp.rTpix(end+1,:)=rTpix;
    tw=4*rTpix+1; th=4*rTpix+1; %Size of support in pixels, derived from si:
    contrast=65; %Contrast
    aspectratio=1; %Aspect ratio
    sc=round(rTpix/2)+1; %Spatial constant of the exponential "hull"
    blobtex=CreateProceduralGaussBlob(wndw,tw,th);
    blobrect=Screen('Rect',blobtex);
    blobpars=[contrast,sc,aspectratio,0];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Gaussian dot-cloud %%%
    sigma=31*[1 0;0 1]; Ndots=65; dtSmooth=30; dtG=.06; countermax=round(30/dtG)+1;
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
        
        sampinitPH(Sexp.sampintervalsec); %pre-allocate memory for sampintervalsec of data
        Sexp.Sdelay(t)= 1; %const delay
        Sexp.Tfixup{t}=GetSecs; %first entry
        pollPH;
        while and(~keyDown,GetSecs-Sexp.Tfixup{t}(end)<Sexp.Sdelay(t)), %Wait for delay (keyDown event breaks loop)
            Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
            Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
            Screen('FillOval',wndw,Sexp.col.RED,Frect+[PHxy(1) PHxy(2) PHxy(1) PHxy(2)]); %show Ftip
            pollPH; Screen('Flip',wndw);
            keyDown=KbCheck;
            
            %%% PRE-CONDITIONING so that Ftip in particular is stationary when target comes up
            othercrit=1;
            if Sexp.jF>Sexp.NthreshPH,
                if std(dist(Sexp.pix2mm(Ftmp(Sexp.jF+[-Sexp.NthreshPH:0],1:2)),EXP(t,1:2),2))>Sexp.FprethreshSD,
                    othercrit=0; end, end
            if or(~othercrit,dist(PHxy(1:2),Sexp.xy0,2)>sz),
                Sexp.jF=0; Sexp.Tfixup{t}(end+1)=GetSecs; else keyDown=1; end %add new start time because finger is not reset properly
            pollPH; end, end %begin reach sequence
    clc; cprintf('blue',['trial # ' num2str(t)  ' \n']);
    
    if ~Sexp.quitExp,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Give countdown to go %%%
        onnow=0; Sexp.t0(t)=GetSecs; Sexp.jF=0; pollPH; %start of recording session
        for n=1:3,
            while GetSecs-Sexp.t0(t)<n/2,
                if and(GetSecs-Sexp.t0(t)>.65,GetSecs-Sexp.t0(t)<1.15), onnow=1;
                    if GetSecs-Sexp.t0(t)>tlistDxy(counter), counter=(counter~=SZdxy3)*counter+1;
                        Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
                        Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
                        Screen('FillOval',wndw,Sexp.col.RED,Frect+[PHxy(1) PHxy(2) PHxy(1) PHxy(2)]); %show Ftip
                        Screen('DrawDots',wndw,Dxy(:,:,counter)',3,Sexp.col.BLU,[],2); %target cloud
                        Screen('Flip',wndw); end
                elseif onnow, onnow=0;
                    Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
                    Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
                    Screen('FillOval',wndw,Sexp.col.RED,Frect+[PHxy(1) PHxy(2) PHxy(1) PHxy(2)]); %show Ftip
                    Screen('Flip',wndw); end
                pollPH; end, end
        Screen('Flip',wndw);
        
        
        %%%%%%%%%%%%%%%%%%%
        %%%%%% REACH %%%%%%
        %%%%%%%%%%%%%%%%%%%
        Sexp.Tnow=EXP(t,3:4);
        goodindex=[0 0 0];
        Sexp.ithreshF=[0 0]; %reset start/end indices for next trial
        startz=PHxy(3); Sexp.Tgo(t)=GetSecs;
        
        playblocking(snd.TINK);
        Sexp.Tbeep(t)=GetSecs; counter=1; tcounter=Sexp.Tbeep(t);
        while Sexp.ithreshF(2)==0, %or(Sexp.ithreshF(2)==0,Sexp.jF<Sexp.ithreshF(2)),
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GET DATA AND UPDATE SCREEN IMAGE %%%
            pollPH; %PHpos=PHxy(1:2)';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% WAITING for the finger mvmt to begin & end %%%
            stdcheck(0); %set Sexp.ithreshF as Ftip mvmt starts and ends
            
            %%%%%%%%%%%%%%%%%%%%%
            %%% Maintain Stim %%%
            if counter>100, counter=1; tcounter=GetSecs; end
            if GetSecs-tcounter>tlistDxy(counter), counter=counter+1;
                if counter==countermax, counter=1; tcounter=GetSecs; end
                Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
                Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
                Screen('DrawDots',wndw,Dxy(:,:,counter)',3,Sexp.col.BLU,[],2); %target cloud
                Screen('FillOval',wndw,Sexp.col.RED,Frect+[PHxy(1) PHxy(2) PHxy(1) PHxy(2)]); %show Ftip
                Screen('Flip',wndw); end; end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% RECORD finger mvmt begin & end data %%%
        if ~BTnow, pollPH
 
        
        
        
%%%% UNCOMMENT IF TARGET SPARKLE PAUSES DURING FEEDBACK BEEP
%         %erase target briefly during FBbeep
%         Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
%         Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
%         Screen('Flip',wndw);
        



        %FB blob
            FBpix=[Ftmp(Sexp.ithreshF(2),1)+Sexp.x0 Sexp.y0-Ftmp(Sexp.ithreshF(2),2)]-Sexp.mm2pix(EXP(t,3:4));
            misshit=IsInCirc(FBpix,[0 0 rTpix]);
            Sexp.Fend=[Sexp.Fend(2:Sexp.NrT,:); FBpix];
            rTpix=round(computeRT(Sexp.Fend,.7,125));
            
            %FB sound
            if misshit, play(snd.GOODBEEP); else play(snd.RASPBERRY); end
            blbcolnow=blbCOL((misshit)+1,:);
            tcounter=GetSecs; counter=1;
            while GetSecs-tcounter<2, pollPH; %show fb for 2s
                if GetSecs-tcounter>tlistDxy(counter), counter=counter+1;
                    Screen('FillOval',wndw,Sexp.col.BLU,startargetrectO); %fixation outer
                    Screen('FillOval',wndw,Sexp.col.BLK,startargetrectI); %fixation inner
                    Screen('DrawTextures',wndw,blobtex,[],CenterRectOnPoint(blobrect,Tpix(1),Tpix(2))',0,[],[],blbcolnow,[],kPsychDontDoRotation,blobpars');
                    Screen('DrawDots',wndw,Dxy(:,:,counter)',3,Sexp.col.BLU,[],2); %end %target cloud
                    Screen('FillOval',wndw,Sexp.col.RED,Frect+[PHxy(1) PHxy(2) PHxy(1) PHxy(2)]); %show Ftip
                    Screen('Flip',wndw); end, end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% UPDATE DATA STORAGE %%%
            PHdat=sampupdatePH(PHdat); %adds a t'th cell to the PHdat cell vector
            
            %%%%%%%%%%%%%%%%%%%%%
            %%% FILL EXP DATA %%%
            %%%%%%%%%%%%%%%%%%%%%
            Sexp.ifinstart(t)=Sexp.ithreshF(1);
            Sexp.ifinend(t)=Sexp.ithreshF(2);
            EXP(t,9:10)=PHdat{t}(Sexp.ifinstart(t),1:2);                      %mvmt strt xy
            EXP(t,11:12)=PHdat{t}(Sexp.ifinend(t),1:2);                       %mvmt end xy
            EXP(t,17:18)=PHdat{t}([Sexp.ifinstart(t) Sexp.ifinend(t)],4); end %mvmt s/e time
        Sexp.BTlist(t)=Sexp.BTlist(t)+BTnow; save tmpPCT Sexp Ftmp PHxy FBpix PHdat;  end, end
Spre=Sexp; datSAVEvg(ELdat,PHdat,EXP); 
