function EyelinkDrawCalTargetTEST(el, x, y)
global Sexp window

% draw simple calibration target
%
% USAGE: rect=EyelinkDrawCalTarget(el, x, y)
%
%		el: eyelink default values
%		x,y: position at which it should be drawn
%		rect: 

% simple, standard eyelink version
%   22-06-06    fwc OSX-ed

%[width, height]=Screen('WindowSize', window);

sz=round(el.calibrationtargetsize/100*Sexp.resx);
inset=round(el.calibrationtargetwidth/100*Sexp.resx);

caltargetrect=CenterRectOnPoint([0 0 sz sz], x, y);
Screen('FillOval', window, Sexp.col.BLU,caltargetrect);
Screen('FillOval',window,Sexp.col.BLK,InsetRect(caltargetrect,inset,inset));
Screen('Flip',window);

% mglClearScreen;
% %mglGluAnnulus(round(x*.61),round(y*.73),2,8,[1 1 1],24);
% mglGluAnnulus(x,y,2,8,[1 1 1],24);
% mglFlush;