function [tlist]=setLissParams(t)
global Sexp
Sexp.Aprime(t,:)=75+20*rand(1,2);
if Sexp.singles==0,
    Sexp.dtSmooth(t)=Sexp.dtSmooth(t-1); Sexp.fXY(t,:)=Sexp.fXY(t-1,:);
    if max(Sexp.lissvar(find(Sexp.lissvar(:,1)>0,1,'last'),:))<Sexp.lissTHRESH,
        if diff(Sexp.fXY(t,:))==0,
            if Sexp.fXY(t,1)<5, Sexp.fXY(t,:)=[5 5]; Sexp.dtSmooth(t)=Sexp.dtSmooth(t-1)+3;
            else Sexp.fXY(t,:)=[7 5]; Sexp.dtSmooth(t)=Sexp.dtSmooth(t-1)+10; end
        else Sexp.dtSmooth(t)=Sexp.dtSmooth(t-1)-2; end
    elseif t>2, Sexp.dtSmooth(t)=Sexp.dtSmooth(t-1)+4; end
else Sexp.dtSmooth(t)=Sexp.dtSmooth(1); Sexp.fXY(t,:)=Sexp.fXY(1,:); end

Nupdate=ceil(Sexp.dtSmooth(t)*Sexp.nominalHertz);
tlist=linspace(0,Sexp.dtSmooth(t),Nupdate);

Sexp.OArot(t)=(diff(Sexp.fXY(t,:))~=0)*pi/8*(2*rand-1);
Sexp.SGN(t,:)=sign(rand(1,2)-.5);
if Sexp.fXY(t,1)==Sexp.fXY(t,2), %circles
    Anow=(1-exp(-tlist(1:end-1)'/.3))*Sexp.Aprime(t,:);
    Sexp.Fpos{t}=rot([[0; Sexp.SGN(end,1)*Anow(:,1).*sin(linspace(0,Sexp.fXY(t,1)*(2*pi),Nupdate-1))'; 0]...
        [0; Sexp.SGN(end,2)*Anow(:,2).*sin(linspace(0,Sexp.fXY(t,2)*(2*pi),Nupdate-1)-pi/2)'; 0]],-Sexp.OArot(t)); %following positions
else %LISS
    Sexp.Fpos{t}=rot([[0; Sexp.SGN(end,1)*Sexp.Aprime(t,1)*sin(linspace(0,Sexp.fXY(t,1)*(2*pi),Nupdate-1))'; 0]...
        [0; Sexp.SGN(end,2)*Sexp.Aprime(t,2)*sin(linspace(0,Sexp.fXY(t,2)*(2*pi),Nupdate-1))'; 0]],-Sexp.OArot(t)); end %following positions
Sexp.Fpos{t}(:,3)=[tlist Sexp.dtSmooth(t)];
%figure(1); plot(Sexp.Fpos{t}(:,1),Sexp.Fpos{t}(:,2),'k.')
