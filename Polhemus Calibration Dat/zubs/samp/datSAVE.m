
function datSAVE(ELdat,PHdat,EXP) %#ok<INUSD,*NASGU> 
global Sexp Ftmp Etmp
if nargin<2, PHdat=[]; end
if nargin<1, ELdat=[]; end
if ~Sexp.quitExp,
    if ~isempty(Ftmp), PHdat=sampupdatePH(PHdat); end %any un-'fixed' PHdata is added to the PHdat cell
    if ~isempty(ELdat), if ~isempty(Etmp{Sexp.EYEnow(1)}), ELdat=sampupdateEL(ELdat); end, end, end %any un-'fixed' ELdata is added to the ELdat cell array

Sexp=rmfield(Sexp,Sexp.FNlist);
if Sexp.Snum<10, Snum=['0' num2str(Sexp.Snum)]; else Snum=num2str(Sexp.Snum); end
if nargin<3, save(['DAT\' Sexp.Stype Snum Sexp.PhaseNow '.mat'],'ELdat','PHdat','Sexp');
else save(['DAT\' Sexp.Stype Snum Sexp.PhaseNow '.mat'],'EXP','ELdat','PHdat','Sexp'); end