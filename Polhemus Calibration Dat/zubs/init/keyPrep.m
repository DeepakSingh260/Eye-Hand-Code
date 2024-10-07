
function keyPrep %#ok<*NASGU,*AGROW>
global snd Sexp

if ismac,
    Sexp.keys.shift=225; Sexp.keys.ctrl=224; Sexp.keys.tab=43; Sexp.keys.capslock=57;
    
    %%% FIND up/down keycodes
    disp('press up arrow'); clear KbCheck; keyCode=0;
    [~,keyCode]=KbStrokeWait; play(snd.CLICK);
    keys.udCode=find(keyCode);
    disp('press down arrow'); clear KbCheck; keyCode=0;
    [~,keyCode]=KbStrokeWait; play(snd.CLICK);
    Sexp.keys.udCode(2)=find(keyCode);
    
    %%% FIND left/right keycodes
    disp('press left arrow'); clear KbCheck; keyCode=0;
    [~,keyCode]=KbStrokeWait; play(snd.CLICK);
    Sexp.keys.lrCode=find(keyCode);
    disp('press right arrow'); clear KbCheck; keyCode=0;
    [~,keyCode]=KbStrokeWait; play(snd.CLICK);
    Sexp.keys.lrCode(2)=find(keyCode);
    
    %%% FIND esc keycode
    disp('press ESC key'); clear KbCheck; keyCode=0;
    [~,keyCode]=KbStrokeWait; play(snd.CLICK);
    Sexp.keys.escCode=find(keyCode);
    
    %%% FIND esc keycode
    Sexp.keys.numCode=nan(1,4);
    for N=1:4, disp(['press ' num2str(N) ' key']); clear KbCheck; keyCode=0;
        [~,keyCode]=KbStrokeWait; play(snd.CLICK);
        Sexp.keys.numCode(N)=find(keyCode); end, end

%PC keycodes
if ispc,  Sexp.keys.shiftCode=16;  Sexp.keys.ctrlCode=17;  Sexp.keys.tabCode=9;  Sexp.keys.capslockCode=20;
    Sexp.keys.escCode=27; Sexp.keys.lrCode=[37 39]; Sexp.keys.udCode=[38 40]; Sexp.keys.numCode=49:57; end

