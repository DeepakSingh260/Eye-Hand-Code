function sampinitPH(Nsecs)
global Ftmp Sexp
if nargin==0, Nsecs=Sexp.sampintervalsec; end

Ftmp=nan(round(1.1*Sexp.PHrate*Nsecs),Sexp.nsampPH);
Sexp.jF=[0];