% logpeakval.m
%
% usage: [logPout]=logpeakval(logPvec);
%
% Computes interpolated peak location given index of sampled peak
%
% INPUTS
% logPvec: [likelihoods abscissaValues]
% optional:
%     ind: index or indices of peak
%     epn: exponent for weighting (bet .25 and .5 works best)
%
% OUTPUT
% peakval: abscissa value corresponding to interpolated peak
%
% teh wrote it. [10.12.11] ** DO NOT DISTRIBUTE

function [peakval ind]=logpeakval(logPvec,ind,epn)
if all(rem(logPvec(:,1),1)==0), logPvec(:,1)=log(logPvec(:,1)/sum(logPvec(:,1))); end
if any(size(logPvec)==1), error('logPvec input must be a 2-column matrix: C1=p or log(p); C2=x'); end
if nargin<3, epn=.375; ind=find(logPvec(:,1)==max(logPvec(:,1))); end
if nargin==2, if all(rem(ind,1)==0), epn=.375; else epn=ind; ind=find(logPvec(:,1)==max(logPvec(:,1))); end, end
if length(ind)==1,
	if or(ind==1,ind==length(logPvec)), peakval=logPvec(ind,2);
	else logPin=logPvec(ind-1:ind+1,1);
		d1=logPin(2)-logPin(1);
		d2=logPin(2)-logPin(3);
		peakval=(d1^epn*logPvec(ind+1,2)+d2^epn*logPvec(ind-1,2))/(d1^epn+d2^epn); end
else peakval=mean(logPvec(ind,2)); end
if nargout==2, ind=round(median(ind)); end