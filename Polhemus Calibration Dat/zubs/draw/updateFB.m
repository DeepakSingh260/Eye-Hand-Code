function updateFB(tnow,h)
global snd Sexp

SoundNow=1;
if tnow>Sexp.Npractice,
    Sexp.dscore=IsInCirc(Sexp.EPmm,[Sexp.Tmm Sexp.rT(2)]);
    drawTarg(Sexp.Tmm,[2*Sexp.rT(2); 2*Sexp.rT(2)-.7], Sexp.FBcol{Sexp.dscore+1},0);  %show Target FB halo
    Sexp.score(h)=Sexp.score(h)+Sexp.dscore;
	Sexp.HR=Sexp.score(h)/(tnow-Sexp.Npractice);
    displayScore(Sexp.HR); end
drawTarg(Sexp.Tmm,Sexp.Tr,Sexp.Tcol,0); %target
drawTarg(Sexp.EPmm,Sexp.EPr,Sexp.EPcol,1); %show result point of last move

if and(tnow>=Sexp.Npractice-1,Sexp.dscore==1), SoundNow=2; end
if SoundNow==1, play(snd.CLICK);
else play(snd.GOODBEEP); end


