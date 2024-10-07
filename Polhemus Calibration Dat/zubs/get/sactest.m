%test for saccade
%return timing and whether it has begun(1) or begun & ended (2)

function sactest
global Etmp Sexp

if Sexp.ithreshE(2,2)==1,
    TFflag=1;
    while TFflag~=0,
        primarysecondary=max([1; sum(Sexp.ithreshE(:,1)>1); sum(Sexp.ithreshE(1,:)>1)]);
        TFflag=stdcheck(1,primarysecondary); %1 is for eye (0 for fin)

        %new saccade event?
        if TFflag, sacstartend=sum(Sexp.ithreshE(primarysecondary,:)>1); %last detect is sac start or end?
            itmp=Sexp.ithreshE(primarysecondary,sacstartend);
            
            %record time and xy-pos of new
            Sexp.sactime(primarysecondary,sacstartend)=Etmp{Sexp.EYEnow(1)}(itmp,4);
            for Enow=Sexp.EYEnow,
                Sexp.sacxy(2*(primarysecondary-1)+sacstartend,:)=Sexp.sacxy(2*(primarysecondary-1)+sacstartend,:)+Etmp{Enow}(itmp,1:2)/Sexp.Neye; end, end, end, end, end %EP
