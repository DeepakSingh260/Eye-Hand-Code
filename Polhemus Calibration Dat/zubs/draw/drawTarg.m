%takes a decreasing set of target radii(in x,y) and a set of corresponding colors
%makes a bullseye with size(rTlist,1) concentric circles
%
%    Tpix: (1x2) vector for target's position on screen (screenCoordinates)
%  rTlist: (Nx2) matrix of diameters (x,y-diams) for filled ovals
%   Clist: (Nx3) matrix of RGB triplets for the filled ovals

function drawTarg(Tin,dlist,Clist,flush)
global wndw Sexp

Tpix=Sexp.mm2pix(Tin(1:2)); %target location
Dpix=Sexp.mm2pix0(dlist); %radii for target
for n=1:size(Dpix,1),
    Screen('FillOval',wndw,Clist(n,:),[Tpix-Dpix(n,:)/2 Tpix+Dpix(n,:)/2]); end
if flush, Screen('Flip', wndw); end