function displayScore(Score)
global wndw Sexp

if ~isempty(Sexp.dscore),
	Screen('TextSize',wndw,20);
	DrawFormattedText(wndw,['Hit Rate: ' num2str(round(100*Score)) '%'],Sexp.resx *.8,Sexp.resy * 0.05,Sexp.col.GRY);
	Screen('TextSize',wndw,Sexp.TextSize); end
