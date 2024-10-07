% rotRel.m
% Finds the signed rotation (rad) of points in the plane relative to other 
% points, such that a rotation of Min by dTheta, Mout=rot(Min,dTheta);, will 
% fall on the radius of the Standard (i.e., irrespective of scaling)
%
% Range of outputs is (-pi,pi): CW/+, CCW/-
%
% usage:
% [dTheta]=rotrel(Min,Standard)
%
% INPUTS
%      Min: Input matrix (eg., data values) arranged as two columns [X Y])
% Standard: Same as above, but for location of standard to which the
%              angle of the Actual will be compared: DEFAULT=[1 0];
% teh wrote it (Oct 2, 2009)

function [dTheta]=rotRel(Min,Standard)
if nargin==1, Standard=[1 0]; end 

%if size of standard or Min is 1, and the other is not, copy it down
if and(size(Min,1)>1,size(Standard,1)==1), Standard=ones(size(Min,1),1)*Standard; end
if and(size(Standard,1)>1,size(Min,1)==1), Min=ones(size(Standard,1),1)*Min; end

ThetaM=atan2(Min(:,1),Min(:,2));
ThetaS=atan2(Standard(:,1),Standard(:,2));
dTheta=ThetaM-ThetaS;
inds=find(and(sign(ThetaM)==-sign(ThetaS),~isnear(abs(dTheta),pi)));
dTheta(isnear(abs(dTheta),pi))=pi;

%When Min and Standard are on opposite sides of vertical meridian of lower hemisphere
indss=inds(and(sign(Min(inds,2))==-1,sign(Standard(inds,2))==-1)); 
if ~isempty(indss), dTheta(indss)=sign(ThetaS(indss)).*(2*pi-sum(abs([ThetaM(indss) ThetaS(indss)]),2)); end

%When Min and Standard are in upperL&lowerR or lowerL&upperR quadrants, resp
indss=inds(~and(sign(Min(inds,2))==-1,sign(Standard(inds,2))==-1));
if ~isempty(indss),
	tmp=[sign(ThetaS(indss)).*(2*pi-sum(abs([ThetaM(indss) ThetaS(indss)]),2)) sign(ThetaM(indss)).*sum(abs([ThetaM(indss) ThetaS(indss)]),2)];
	for n=1:length(indss), dTheta(indss(n))=tmp(n,abs(tmp(n,:))==min(abs(tmp(n,:)),[],2)); end, end, end
