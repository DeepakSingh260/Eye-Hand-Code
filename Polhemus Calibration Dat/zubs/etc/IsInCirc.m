% inside=IsInCirc([x y],[x y rad])
%
% Decide if location x,y is inside of the passed rect.
%
% 3/5/97  teh  Wrote it.
function inside=IsInCirc(coords,Circ)
if size(Circ,2)==4, Circ=[mean(Circ(:,[1 3]),2) mean(Circ(:,[2 4]),2) diff(Circ(:,[1 3]),1,2)./2]; end
inside=sqrt(sum((coords-Circ(:,1:2)).^2,2))<Circ(:,3);
