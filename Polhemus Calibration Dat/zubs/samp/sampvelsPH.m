function vout=sampvelsPH
global Ftmp Sexp

vout=sum(Ftmp(1:Sexp.jF,5:7).^2);