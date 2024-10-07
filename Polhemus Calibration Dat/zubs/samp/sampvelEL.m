

function [velnow]=sampvelEL
global Etmp Sexp
ELtmp=zeros(2,3);
for Enow=Sexp.EYEnow, 
    ELtmp=ELtmp+Etmp{Enow}(Sexp.jE+[-1 0],[1 2 4]); end
velnow=sampfixEL(ELtmp/Sexp.Neye);