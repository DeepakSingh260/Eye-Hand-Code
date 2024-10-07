%test for fmovecade
%return timing and whether it has begun fmovedone=ind1, or begun & ended fmovedone=[ind1,ind2]

function [fmovedone]=fmovetest(fmovedone) %for Ftip, have fmovedone be the time index of condition true
global Ftmp Sexp
%fmovetime=nan(1,2);

v=abs(diff(dist(Sexp.Fnow,Sexp.pix2mm(Ftmp(1:Sexp.jF,1:2)),Sexp.Tnow-Sexp.Fnow)))./diff(Ftmp(1:Sexp.jF,4));
switch Sexp.NvthreshPH(length(fmovedone)+1), %allow for easy switching among strictnesses of threshold criterion
	case 2, vlist=[v(1:end-1) v(2:end)];
	case 3, vlist=[v(1:end-2) v(2:end-1) v(3:end)];
	case 4, vlist=[v(1:end-3) v(2:end-2) v(3:end-1) v(4:end)];
	case 5, vlist=[v(1:end-4) v(2:end-3) v(3:end-2) v(4:end-1) v(5:end)];
	case 6, vlist=[v(1:end-5) v(2:end-4) v(3:end-3) v(4:end-2) v(5:end-1) v(6:end)];
	case 7, vlist=[v(1:end-6) v(2:end-5) v(3:end-4) v(4:end-3) v(5:end-2) v(6:end-1) v(7:end)];
	case 8, vlist=[v(1:end-7) v(2:end-6) v(3:end-5) v(4:end-4) v(5:end-3) v(6:end-2) v(7:end-1) v(8:end)]; end

%have to pre-condition the post-go data so that it starts with stationary finger
%with combined var criterion and vel criterion

switch length(fmovedone),
	case 0, 
		fmovedone=find(all(vlist>Sexp.PHvthresh,2),1); %has fmove started?
		Sexp.E15=Ftmp(fmovedone,4);
	case 1,
		vlist(1:fmovedone)=nan;
		iend=find(all(vlist<Sexp.PHvthresh,2),1); %has fmove ended?
		if ~isempty(iend), fmovedone(2)=iend; end
	case 2, end %bypass

