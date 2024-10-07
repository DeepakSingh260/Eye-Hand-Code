%sampupdate adds the Nth row to the cell array sampnow (which is an Nx2 array)
function [sampnew]=sampupdatePH(sampold)
global Ftmp Sexp
Ftmp=[Ftmp(:,1)+Sexp.x0 Sexp.y0-Ftmp(:,2) Ftmp(:,3:end)];
if nargin==0, sampnew=cell(1); 
else sampnew=[sampold; cell(1)]; end

sampnew{end}=sampfixPH(Ftmp(1:Sexp.jF,:));
Ftmp=[];