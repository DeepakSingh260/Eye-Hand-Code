%stdcheck check std of eye and hand

function tflogical=stdcheck(FEflag,primarysecondary)
global Ftmp Etmp Sexp

tflogical=0;
if and(~FEflag,Sexp.jF>Sexp.NthreshPH), %engage fingertip check
    startend=sum(Sexp.ithreshF>0)+1;
    if and(Sexp.jF>Sexp.NthreshPH,startend<3),
        Fdat=Sexp.pix2mm(Ftmp(Sexp.jF+[-Sexp.NthreshPH:0],1:2));
        if Sexp.ithreshF(1)==0, tflogical=std(dist(Fdat,Sexp.Fnow,2))>Sexp.FthreshSD;
        else tflogical=std(dist(Fdat,Fdat(end,:),2))<Sexp.FthreshSD*2; %std of ftip dist from fixation target (Fnow) less than Sexp.FthreshSD mm
            if tflogical, [keyDown,~,kbNow]=KbCheck;
                if keyDown, tflogical=find(kbNow)~=Sexp.keys.shiftCode; end, end
            if tflogical,
                tflogical=dist(Fdat(end,:),Sexp.Fnow,2)>Sexp.ReachDist/3; end, end
        if tflogical, Sexp.ithreshF(startend)=Sexp.jF-Sexp.NbackPH; end, end %this test occurs online, eyetest occurs offline (or in chuncks)
    
elseif and(FEflag,Sexp.jE>Sexp.NthreshEL+1), %engage eye check
    if and(Sexp.ithreshE(2,2)==1,Sexp.jE>Sexp.NthreshEL),
        Edat=Sexp.pix2mm(Etmp{Sexp.EYEnow(1)}(1:Sexp.jE,1:2))/Sexp.Neye;
        if Sexp.Neye==2, Edat=Edat+Sexp.pix2mm(Etmp{2}(1:Sexp.jE,1:2))/2; end
        Edat=dist(Edat,2);
        sacstartend=sum(Sexp.ithreshE(primarysecondary,:)>1)+1;
        
        LEdat=max([Sexp.ithreshE(:)'+Sexp.NbackEL(:)' find(Etmp{Sexp.EYEnow(1)}(:,4)<Sexp.ELgo(end),1,'last')]); dnow=nan(Sexp.jE,1);
        LEdat=max([LEdat Sexp.NthreshEL+1]);
        for n=LEdat:Sexp.jE, dnow(n)=std(Edat(n+[-Sexp.NthreshEL:0])); end %sliding window of std vals
        if sacstartend==1, %if start unset (no index currently set into ithreshE,1)
            tflogical=find(dnow>Sexp.EthreshSD(primarysecondary)); %start occurs when sd first rises above thresh for two timesteps
        else %if end unset (ithreshE,2 is still default val)
            tflogical=find(dnow<Sexp.EthreshSD(primarysecondary)); end %end occurs when sd first falls below thresh for two timesteps
        if length(tflogical)>1, tflogical=tflogical(find(diff(tflogical)==1,1,'first')); else tflogical=[]; end
        %if start found, set into ithreshE,1 (minus NbackEL) and set flag to 1
        if ~isempty(tflogical), Sexp.ithreshE(primarysecondary,sacstartend)=max([1+sacstartend tflogical-Sexp.NbackEL(primarysecondary)]); tflogical=1;
        else tflogical=0; end, end, end %otherwise set flag to 0
