%test for saccade
%return timing and whether it has begun(1) or begun & ended (2)

function [sacstartend sactime sacxy]=sactest(primarysecondary)
global Etmp Sexp
sacstartend=0; sactime=zeros(1,2); sacxy=zeros(2);

v=nan(Sexp.jE-1,2);
for Enow=Sexp.EYEnow,
	sampout=Sexp.pix2mm(Etmp{Enow}(1:Sexp.jE,1:2));
	if primarysecondary==1,
		v(:,Enow)=windowmean(abs(diff(dist(Sexp.Fnow,sampout,Sexp.Tnow-Sexp.Fnow))),31,[])./diff(Etmp{Enow}(1:Sexp.jE,4));
	else v(:,Enow)=windowmean(abs(diff(dist(sampout,Sexp.Fnow,2))),31,[])./diff(Etmp{Enow}(1:Sexp.jE,4)); end, end

v=nanmean(v,2);
switch Sexp.NvthreshEL(primarysecondary), %allow for easy switching among strictnesses of threshold criterion
	case 2, vlist=[v(1:end-1) v(2:end)];
	case 3, vlist=[v(1:end-2) v(2:end-1) v(3:end)];
	case 4, vlist=[v(1:end-3) v(2:end-2) v(3:end-1) v(4:end)];
	case 5, vlist=[v(1:end-4) v(2:end-3) v(3:end-2) v(4:end-1) v(5:end)];
	case 6, vlist=[v(1:end-5) v(2:end-4) v(3:end-3) v(4:end-2) v(5:end-1) v(6:end)];
	case 7, vlist=[v(1:end-6) v(2:end-5) v(3:end-4) v(4:end-3) v(5:end-2) v(6:end-1) v(7:end)];
	case 8, vlist=[v(1:end-7) v(2:end-6) v(3:end-5) v(4:end-4) v(5:end-3) v(6:end-2) v(7:end-1) v(8:end)]; end
vlist(1:Sexp.isacend,1)=nan;
vlist(Etmp{Enow}(:,4)<Sexp.ELgo(end),1)=nan;

inow=find(all(vlist>Sexp.ELvthresh(primarysecondary),2),1); %has saccade started?
if ~isempty(inow), sacstartend=1;
	sactime(1)=Etmp{Sexp.EYEnow(1)}(inow,4);
	vlist(1:inow,1)=nan;
	iend=find(all(vlist<Sexp.ELvthresh(2),2),1); %has saccade ended?
	if ~isempty(iend), inow(2)=iend; Sexp.isacend=iend;
		sactime(2)=Etmp{Sexp.EYEnow(1)}(iend,4); sacstartend=2; end
	
	for n=1:sacstartend,
		for Enow=Sexp.EYEnow,
			sacxy(n,:)=sacxy(n,:)+Etmp{Enow}(inow(n),1:2); end
		sacxy(n,:)=sacxy(n,:)/Sexp.Neye; end, end %EP
