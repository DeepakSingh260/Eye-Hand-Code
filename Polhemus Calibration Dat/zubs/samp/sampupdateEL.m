%sampupdate adds the Nth row to the cell array Etmp (which is an Nx2 array)
function [sampnew]=sampupdateEL(sampold)
global Etmp
sampnew=sampfixEL(Etmp);

if nargin==1, sampnew=[sampold; sampnew]; end
Etmp=cell(1,2);