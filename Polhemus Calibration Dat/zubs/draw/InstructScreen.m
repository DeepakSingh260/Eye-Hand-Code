function InstructScreen(instructcell,centerrelheight,colnow)
global wndw
if nargin<2, error('Must provide 2 inputs'); end
if length(instructcell)~=length(centerrelheight), error('''instructcell'' and ''centerrelheight'' inputs do not match in size'); end

for n=1:length(instructcell),
DrawFormattedText(wndw, instructcell{n}, 'center', centerrelheight(n), colnow{n}); end
Screen('Flip', wndw); KbStrokeWait; %wait for buttonpress to continue

