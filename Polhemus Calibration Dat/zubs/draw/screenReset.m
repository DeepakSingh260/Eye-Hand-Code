function screenReset
global wndw %Sexp

for n=1:2,
	%Screen('FillRect',wndw,Sexp.BGcolor); % Clear buffer with black.
	Screen('Flip',wndw); end % blank screen.

