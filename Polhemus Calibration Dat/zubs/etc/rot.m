% rot.m
% usage:
% [Mout]=rot(Min,theta,IfDeg)
%
% Min: Input matrix (eg., data values) arranged as two columns)
% theta: rotation angle (radians; otherwise set IfDeg flag)
% rot matrix: [cos(theta) -sin(theta)
%              sin(theta)  cos(theta)]
%
% The rotation direction uses the atan2 convention: +/CCW, -/CW
%
% NOTE: There is a massive speed improvement when all rows of Min are
% rotated by the same angle.
% teh wrote it [Oct 2, 2009]

function [Mout]=rot(Min,theta,IfDeg,ScaleFactor) %#ok<*AGROW>
if nargin<3, IfDeg=0; end
if IfDeg, theta=pi*theta./180; end
if nargin<4, ScaleFactor=1; end

if all(theta==0), 
	if length(ScaleFactor)==1, ScaleFactor=ScaleFactor*ones(size(Min,1),1); end
	Mout=Min.*(ScaleFactor(:)*[1 1]);
else     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MAIN LOOP (once for each dim3) %%%
    Morig=Min;
    for S=1:size(Morig,3), Min=Morig(:,:,S);
        % Check that matrices were formatted properly
        if and(size(Min,1)==2,size(Min,2)~=2), Min=Min';
        elseif all(size(Min)~=2), error('Input matrix should be Nx2'); end
		if length(ScaleFactor)==1, ScaleFactor=ScaleFactor*ones(size(Min,1),1); end
		Min=Min.*(ScaleFactor(:)*[1 1]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Single rotation of all points %%%
        if length(theta)==1,
            Mout(:,:,S)=Min*[cos(theta) sin(theta); -sin(theta)  cos(theta)];
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Multiple rotations of a single point %%%            
            if size(Min,1)==1, Mout(:,:,S)=(ones(length(theta),1)*Min); inds=1:length(theta);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Differrent rotations for each point %%%
            else Mout=Min;
                inds=find(theta~=0); end
            for n=1:length(inds), %only rotate if theta~=0
                Mout(inds(n),:,S)=Mout(inds(n),:)*[cos(theta(inds(n))) sin(theta(inds(n))); -sin(theta(inds(n)))  cos(theta(inds(n)))]; end, end, end, end