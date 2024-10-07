function SN=xtractSnum(CP)
CPstr='CP';
if ~exist('DAT','dir'), mkdir('DAT'); SN=1; else
a=what('DAT'); %extract next Snum of C or P type
   SN=1; SNstr='01';
    while sum(strcmpi([CPstr(CP+1) SNstr 'dat.mat'],a.mat))==1,
        SN=SN+1; if SN<10, SNstr=['0', num2str(SN)]; else SNstr=num2str(SN); end, end, end