function retrievetempsubdata
global Sexp
LRstr={'L','R'}; CPstr={'C','P'}; LMstr={'L','M'};
load tmpsubdata.mat
Sexp.Snum=options.Snum;
Sexp.Stype=CPstr{options.CP+1};
Sexp.LRhand=LRstr{options.LR+1};
Sexp.LRdmg=''; Sexp.LMa=''; Sexp.FMS=0;
if strcmp(Sexp.Stype,'P'),
    Sexp.LRdmg=LRstr{options.LRdmg+1};
    Sexp.LMa=LMstr{options.LMa+1};
    Sexp.FMS=options.FMS; end
delete tmpsubdata.mat