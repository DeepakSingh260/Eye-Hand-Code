function initFtmp(Dlength, Ncarry)
global Ftmp Sexp 
if nargin==1, Ncarry=0; end

Ftmp=[Ftmp(Sexp.jF-Ncarry+1:Sexp.jF,:); nan(Dlength,7)];
Sexp.jF=Ncarry;
