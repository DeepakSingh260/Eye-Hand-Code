function doELdc(el,t)
global Sexp
clc
cprintf('\n'); cprintf('\n'); 
cprintf('*RED','''[ESC]'''); cprintf('blue',': to accept eye position for drift correction \n');

clear KbCheck
Eyelink('StopRecording');
Eyelink('DriftCorrStart',Sexp.x0,Sexp.y0,0,1,0);
EyelinkClearCalDisplay(el);
[~, DCmessage]=Eyelink('CalMessage'); Sexp.DCdat{t,1}=DCmessage;
Eyelink('ApplyDriftCorr');
Eyelink('StartRecording');
clc
