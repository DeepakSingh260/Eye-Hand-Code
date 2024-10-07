function sd=getFsd
global Sexp Ftmp


sd=std(dist(Sexp.pix2mm(Ftmp(Sexp.jF+[-Sexp.NthreshPH:0],1:2)),2));