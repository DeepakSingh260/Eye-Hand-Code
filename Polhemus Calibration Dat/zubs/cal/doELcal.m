function doELcal(el)
global Sexp
clc
cprintf('\n'); cprintf('\n'); cprintf('*RED','''  c'''); cprintf('blue',': to calibrate\n');
cprintf('*RED','''  v'''); cprintf('blue',': to validate calibration\n');
cprintf('*RED','''[esc]'''); cprintf('blue',': to return to setup '); cprintf('*RED','(2x to start exp)\n');

clear KbCheck
Eyelink('StopRecording');
Eyelink('NewestFloatSample');
EyelinkDoTrackerSetup(el);
Eyelink('StartRecording');    
while Eyelink('NewFloatSampleAvailable')==0, end
evt=Eyelink('NewestFloatSample');
Sexp.EYEnow=[];
while isempty(Sexp.EYEnow),
    Sexp.EYEnow=find(and(abs(evt.pa)>0,abs(evt.pa)<32768)); end %1:Leye, 2:Reye
for Enow=1:2, Sexp.Ccoords{1,Enow}=nan(1,2); end
for Enow=Sexp.EYEnow, Sexp.Ccoords{1,Enow}=[0 0]; end
Sexp.Neye=length(Sexp.EYEnow); clc