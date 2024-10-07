% pollPH.m
%
% usage: pollPH;
% ringbuffer: ring buffer contents [C1:X, C2:Y, C3:Z]
%          T: transformation from polhemus-native to pixel coordinates
%       Ftmp: temp data matrix [x y z t vx vy vz]
% 
% teh wrote it. 3.14.15

function pollPH
global ringbuffer Ftmp PHxy Sexp Tx
newSample=10*double(cell2mat(poldata('127.0.0.1',7234)));             %extract latest raw sample (mm)

if Sexp.PHnumsampavg>1,
    if or(Sexp.jF==0,isempty(ringbuffer)),
        ringbuffer=nan(2,3);                                              %create ring
        Sexp.jF=1; Ftmp=nan(1,Sexp.nsampPH);                              %j index = 1 (create Ftmp)
        ringbuffer(end,:)=Tx(newSample(Sexp.Clist));                      %Sexp.Clist are polhemus native XYZ data columns
        Ftmp(Sexp.jF,1:4)=[nanmedian(ringbuffer,1)-[0 0 Sexp.z0] GetSecs];%add to Ftmp
        %if ~Sexp.corrok, ringbuffer(end,:)=nan; end                      %reset datastream after each Ftmp entry if ~corrok
        
    else %separated for tiny speed boost
        Sexp.jF=Sexp.jF+1; Ftmp(Sexp.jF,:)=nan;                                                         %update j index
        ringbuffer(1:end-1,:)=ringbuffer(2:end,:);                                                      %shift ring
        ringbuffer(end,:)=Tx(newSample(Sexp.Clist));                                                    %Sexp.Clist are polhemus native XYZ data columns
        Ftmp(Sexp.jF,:)=[nanmedian(ringbuffer,1)-[0 0 Sexp.z0] GetSecs zeros(1,3)]; end                 %add to Ftmp (screen-center-origin pixels, up,right-positive)
        %if ~Sexp.corrok, ringbuffer(end,:)=nan; end, end                                               %reset datastream after each Ftmp entry if ~corrok
        %Ftmp(Sexp.jF,5:7)=diff(Ftmp(Sexp.jF-1:Sexp.jF,1:3),1,1)/diff(Ftmp(Sexp.jF-1:Sexp.jF,4),1,1); end %3dims of speed vals
else
    if Sexp.jF==0, Sexp.jF=1;                            %j index = 1 (create Ftmp)
        ringbuffer=Tx(newSample(Sexp.Clist));            %Sexp.Clist are polhemus native XYZ data columns
        Ftmp(1,1:4)=[ringbuffer-[0 0 Sexp.z0] GetSecs];  %add to Ftmp
        
    else %separated for tiny speed boost
        Sexp.jF=Sexp.jF+1; Ftmp(Sexp.jF,:)=nan;                                                               %update j index
        ringbuffer=Tx(newSample(Sexp.Clist));                                                                 %Sexp.Clist are polhemus native XYZ data columns
        Ftmp(Sexp.jF,1:4)=[ringbuffer-[0 0 Sexp.z0] GetSecs]; end, end                                        %add to Ftmp (screen-center-origin pixels, up,right-positive)
        %Ftmp(Sexp.jF,5:7)=diff(Ftmp(Sexp.jF-1:Sexp.jF,1:3),1,1)/diff(Ftmp(Sexp.jF-1:Sexp.jF,4),1,1); end, end %3dims of speed vals
PHxy=[Ftmp(Sexp.jF,1)+Sexp.x0 Sexp.y0-Ftmp(Sexp.jF,2) Ftmp(Sexp.jF,3:end)];

