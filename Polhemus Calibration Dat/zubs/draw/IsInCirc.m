% inside=IsInCirc([x y],[x y R])
%       or ([x y],[x1 y1 x2 y2])
%
% Decide if location x,y is inside of the passed radius.
%
% 3/5/97  teh  Wrote it.
function inside=IsInCirc(coords,Circ)

if size(Circ,2)==4,
    %Circ=[mean(Circ(:,[1 3]),2) mean(Circ(:,[2 4]),2) diff(Circ(:,[1 3]),1,2)./2];
    
    inside=and(and(coords(1)>Circ(1),coords(1)<Circ(3)),and(coords(2)>Circ(2),coords(2)<Circ(4)));
else
    if size(coords,2)~=2, error('Arrange data in two columns'); end
    if and(size(coords,1)>1,size(Circ,1)==1),
        Circ=ones(size(coords,1),1)*Circ; end
    if size(Circ,2)==1, Circ=[zeros(size(Circ,1),2) Circ]; end
    if and(size(coords,1)==1,size(Circ,1)>1), coords=ones(size(Circ,1),1)*coords; end
    
    inside=sqrt(sum((coords-Circ(:,1:2)).^2,2))<Circ(:,3);
end