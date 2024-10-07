function sampinitEL(Nsecs)
global Etmp Sexp
if nargin==0, Nsecs=Sexp.sampintervalsec; end

Etmp=cell(1,2);
if any(Sexp.EYEnow==1), Etmp{1}=nan(round(1.1*Sexp.ELrate*Nsecs),Sexp.nsampEL); end %ELrate*Nsecs samples x 4 cols (xpix,ypix,pdiam,time)
if any(Sexp.EYEnow==2), Etmp{2}=nan(round(1.1*Sexp.ELrate*Nsecs),Sexp.nsampEL); end
Sexp.jE=0;

initELold; 