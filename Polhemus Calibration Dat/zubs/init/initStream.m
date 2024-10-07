function Scal=initStream(Pcell,sampmove,Scal) %#ok<*NASGU> %finds the columns corresponding to x-y plane. assumes the native coordinates are aligned to this plane.
global snd

%eliminate all dimensions that do not at some point satisfy abs(x)>pi
counter=0; goodtogo=0; while ~goodtogo, disp('try-loop')
    try goodtogo=1; counter=counter+1;
        if counter<11,
            disp(['counter: ' num2str(counter)])
            sampmove=cell2mat(sampmove(:));
            mvmax=max(abs(sampmove(:,1:end-1)-ones(size(sampmove,1),1)*nanmedian(sampmove(:,1:end-1))));
            cartdims=find(mvmax>1.05*pi); %allow for 5% error [ignore any bad channels with var=0]
            if length(cartdims)<3, Scal.calgood=0;
                clc; cprintf('*RED',['program has identified ' num2str(length(cartdims)) ' Cartesian dimensions: repeating table calibration \n']); play(snd.RASPBERRY);
            else Scal.calgood=1;
                
                %pairs of x- and y-axis-parallel vectors
                xlist=[2 1; 3 2; 5 4; 6 5; 8 7; 9 8];
                ylist=[1 4; 4 7; 2 5; 5 8; 3 6; 6 9];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% FIND XY PLANE INPUTS %%%
                %Scal.list=nchoosek(find(Scal.vars>0),2);
                Scal.list=nchoosek(cartdims,2);
                Scal.Nxs=nan(size(xlist,1),size(Pcell{1},2)); Scal.Nys=Scal.Nxs; Nx=nan(size(Scal.list,1),2); Ny=Nx; DV=nan(size(Scal.list,1),3);
                for nnx=1:size(xlist,1),
                    Scal.Nxs(nnx,:)=nanmedian(Pcell{xlist(nnx,1)},1)-nanmedian(Pcell{xlist(nnx,2)},1); end
                for nny=1:size(ylist,1),
                    Scal.Nys(nny,:)=nanmedian(Pcell{ylist(nny,1)},1)-nanmedian(Pcell{ylist(nny,2)},1); end
                
                for r=1:size(Scal.list,1),
                    Nx(r,:)=median(Scal.Nxs(:,Scal.list(r,:)),1);
                    Ny(r,:)=median(Scal.Nys(:,Scal.list(r,:)),1);
                    if ~or(any(all(Scal.Nxs(:,Scal.list(r,:))==0,2)),any(all(Scal.Nys(:,Scal.list(r,:))==0,2))),
                        DV(r,:)=[std(rotRel(Scal.Nxs(:,Scal.list(r,:)),Nx(r,:)),1) std(rotRel(Scal.Nys(:,Scal.list(r,:)),Ny(r,:)),1) abs(pi/2-abs(rotRel(Nx(r,:),Ny(r,:))))]; end, end %last term is mean-sq-error
                Scal.imin=find(isnear(dist(DV,2),min(dist(DV(dist(DV,2)>0,:),2)))); %can't be ==0
                Scal.Clist=Scal.list(Scal.imin,:);
                Scal.T=[[Nx(Scal.imin,:) 0]/dist(Nx(Scal.imin,:),2); [Ny(Scal.imin,:) 0]/dist(Ny(Scal.imin,:),2); [0 0 1]]; %rotation tensor
                
                %%%%%%%%%%%%%%%%%%%%%%%
                %%% ADD Z-DIM INPUT %%%
                for n=1:2, irem(n)=find(cartdims==Scal.Clist(n)); end
                cleft=RFL(cartdims,irem);
                
                %all done if length(cleft)==1. otherwise, remove imposters:
                if length(cleft)>1, test=zeros(1,length(cleft));
                    for n=1:length(cleft), test(n)=or(all(sampmove(:,cleft(n))>0),all(sampmove(:,cleft(n))<0)); end
                    if or(sum(test)==0,sum(test>1)), for n=1:length(test), test(n)=and(mvmax(cleft(n))>0,mvmax(cleft(n))<1000); end, end %prefer not to use this criterion, but will do in a pinch
                    cleft=RFL(cleft,find(~test));
                    if length(cleft)~=1, play(snd.RASPBERRY);
                        for n=1:3,
                            clc; pause(.05);
                            cprintf('*RED','Too many raw data columns consistent with the Cartesian dimensions: repeating table calibration \n \n');
                            cprintf('*BLUE','[check that the cube is flush with the underside of the table \n'); pause(.5); end
                        Scal.calgood=0; end, end
                Scal.Clist(3)=cleft(1);  %z-dim is the remaining item in list (need the 1 to prevent crash if multiple matches require repeat of tablecal
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% finalize T and store  center and button positions %%%
                Scal.T(end,end)=sign(median(sampmove(:,cleft(1))));
                Scal.F=(Scal.T*nanmedian(Pcell{end}(:,Scal.Clist))')'; %rotated coordinates into Cartesian tabletop plane (ftip button)
                for n=1:length(Pcell)-1, Scal.P5(n,:)=(Scal.T*nanmedian(Pcell{n}(:,Scal.Clist))')'; end
                for nnx=1:size(xlist,1),
                    Scal.Dxy(nnx,1)=abs(Scal.P5(xlist(nnx,1),1)-Scal.P5(xlist(nnx,2),1));
                    Scal.Dxy(nnx,2)=abs(Scal.P5(ylist(nnx,1),2)-Scal.P5(ylist(nnx,2),2)); end
                Scal.C=mean(Scal.P5,1); end %(Scal.T*mean(Scal.P5,1)')'; end %rotated coordinates into Cartesian tabletop plane ('screen' center)
        else play(snd.RASPBERRY);
            for n=1:3,
                clc; pause(.05);
                cprintf('*RED','Data cannot be made into a regular grid: repeating table calibration \n \n');
                cprintf('*BLUE','[check that the cube is flush with the underside of the table \n'); pause(.5); end
            Scal.calgood=0; end
    catch, goodtogo=0; end, end
