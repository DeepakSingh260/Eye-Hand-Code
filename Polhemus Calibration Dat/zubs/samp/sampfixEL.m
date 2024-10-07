function [sampout]=sampfixEL(sampin)
global Sexp ELxy
dostuff=1;
%if there are no inputs, or if the inputs are integers (i.e., eye-numbers), do CASE 1
if nargin==0, dimlist=1:2; end
if iscell(sampin), dostuff=2;
elseif all(size(sampin)==[2 3]), dostuff=3;
else dimlist=sampin; end %previously had internal clause: if all(rem(sampold,1)==0), dostuff=1; end, but this seemed time-consuming for low-freq. event
switch dostuff
	%sampfix is called with no inputs, or with integer list of dims to send as output
	case 1, sampout=Sexp.pix2mm(ELxy(1:2));
		sampout=sampout(dimlist); %Onow is the output dimension (x:1, y:2); note sampnew has #elements = length(sampold)

    %sampfix is called with the Etmp global as sampin input	
	case 2, sampout=cell(size(sampin)); %sets the size for new dataset
		for Enow=Sexp.EYEnow, %add a column for velocities (in target direction, or in mvmt direction for Sac & Liss, resp)
            sampout{Enow}=sampin{Enow}(1:Sexp.jE,:);
			sampout{Enow}(:,1:2)=Sexp.pix2mm(sampout{Enow}(:,1:2));
			%sampout{Enow}(:,3)=sampout{Enow}(:,5); %pupil area is unchange
            %vels (sampout{}(:,5) not yet filled (2D see below)
            %times (sampout{}(:,4) also unchanged
            if isfield(Sexp,'Tnow'),
                sampout{Enow}(2:Sexp.jE,5)=diff(dist(Sexp.Fnow,sampout{Enow}(:,1:2),Sexp.Tnow-Sexp.Fnow))./diff(sampout{Enow}(:,4)); 
            else sampout{Enow}(2:Sexp.jE,5)=abs(diff(dist(sampout{Enow}(:,1:2),Sexp.Fnow,2)))./diff(sampout{Enow}(:,4)); end, end
	
    %sampfix is called with 2x3 input (current & prev samps [ELxy; ELold]
	case 3, sampin(:,1:2)=Sexp.pix2mm(sampin(:,1:2)); %sampin(:,3)=sampin(:,3);
		if isfield(Sexp,'Tnow'),
			sampout=abs(diff(dist(sampin(:,1:2),Sexp.Fnow,Sexp.Tnow-Sexp.Fnow)))./diff(sampin(:,3)); %in the direction of the target
		else sampout=abs(diff(dist(sampin(:,1:2),Sexp.Fnow,2)))./diff(sampin(:,3)); end, end         %in the direction of the mmvmt
