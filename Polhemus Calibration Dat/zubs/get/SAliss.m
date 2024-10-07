%SAliss

load C1LFdat
t=6; lisstest(ELdat(t*2,:),PHdat{t*2},t);
load tmpdat

figure(1); clf
for s=1:2,
    subplot(2,1,s); hold on
    plot(tF,Fpos(:,s),'k.');
    plot(tF,-Fpos(:,s),'.','Color',.6*[1 1 1]); end
subplot(2,1,1); plot(tF,tmpxE,'b.'); plot(tF,-tmpxF,'c.'); r=axis; axis([0 tF(end) r(3:4)])
subplot(2,1,2); plot(tF,tmpxE,'b.'); plot(tF,-tmpxF,'c.'); r=axis; axis([0 tF(end) r(3:4)])

figure(2); clf
subplot(2,1,1)
[height1 pos1]=histnew((dyE),31,1,'k',[.1 .4 .1]);
hold on; [height2 pos2]=histnew((dxE),31,1,'k',[.25 .6 .25]);
axis([min([pos1(:); pos2(:)]) max([pos1(:); pos2(:)]) 0 max([height1(:); height2(:)])])

subplot(2,1,2)
[height1 pos1]=histnew((dyF),31,1,'k',[.1 .4 .1]);
hold on; [height2 pos2]=histnew((dxF),31,1,'k',[.25 .6 .25]);
axis([min([pos1(:); pos2(:)]) max([pos1(:); pos2(:)]) 0 max([height1(:); height2(:)])])

figure(3); clf
subplot(2,1,1); hold on; 
[hs xs]=histnew(dist([dxE(:) dyE(:)],2)/25,31,1,'k',[.1 .2 .4]); r=axis; 
stat=min([Sexp.lissvar(t,1) xs(hs==max(hs))]); title(num2str(stat))
plot(stat*[1 1],r(3:4),'c--','LineWidth',2); 

subplot(2,1,2); hold on; 
[hs xs]=histnew(dist([dxF(:) dyF(:)],2)/25,31,1,'k',[.1 .2 .4]); r=axis; 
stat=min([Sexp.lissvar(t,2) xs(hs==max(hs))]); title(num2str(stat))
plot(stat*[1 1],r(3:4),'c--','LineWidth',2); 

figure(4); clf
plotmat(Fpos(:,1:2),'k.')