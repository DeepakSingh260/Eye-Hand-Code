function [sampout]=sampfixPH(dimout,dims2)
global Sexp PHxy
if nargin==0, dimout=1:2; sampout=Sexp.pix2mm(PHxy(1:2)); %no input returns xy mm
elseif all(size(dimout)==[2 4]), %2 x 4 input is [PHxy; PHxyOLD], and returns vel in mm/s
    sampout=[Sexp.pix2mm(dimout(:,1:2)) dimout(:,3:4)]; dimout=1;
    if isfield(Sexp,'Tnow'), sampout=[diff(dist(Sexp.Fnow,sampout(:,1:2),Sexp.Tnow-Sexp.Fnow))./diff(sampout(:,4))];
		else sampout(:,end+1)=[diff(dist(Sexp.Fnow,sampout(:,1:2),2))./diff(sampout(:,4))]; end
elseif size(dimout,1)>1, %input is Ftmp matrix
	sampout=dimout(1:Sexp.jF,:); sampout(:,1:2)=Sexp.pix2mm(sampout(:,1:2));
	if nargin==2, dimout=dims2; else dimout=1:size(dimout,2); end
	if max(dimout)>size(sampout,2),
		if isfield(Sexp,'Tnow'), sampout(:,end+1)=[nan; diff(dist(Sexp.Fnow,sampout(:,1:2),Sexp.Tnow-Sexp.Fnow))./diff(sampout(:,4))];
		else sampout(:,end+1)=[nan; diff(dist(Sexp.Fnow,sampout(:,1:2),2))./diff(sampout(:,4))]; end, end, end
sampout=sampout(:,dimout);