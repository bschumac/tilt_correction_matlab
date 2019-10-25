function [U1,V1,W1,varargout] = tiltCorrection(u,v,w,varargin)
% [U1,V1,W1] = tiltCorrection(u,v,w,varargin) corrects the wind
% velocity components from non-zero tilt angles and computes them in the wind-based
% coordinate system. The tilt angle correction is done using the Planar Fit
% (PF) method [1].
%
% Notations:
% M: number of samples
% N: Number of time step
%
% Input:
% u: [M x N] matrix: First wind component recorded by a sonic anemometer (First horizontal component)
% v: [M x N] matrix: Second wind component recorded by a sonic anemometer (Second horizontal component)
% w: [M x N] matrix: Third wind component recorded by a sonic anemometer (vertical component)
% varargin:
% method: 'PF' to use [1].
%
% Output:
% U1: [M x N] matrix: Along wind component
% V1: [M x N] matrix: Across wind component
% W1: [M x N] matrix: Vertical wind component
% varargout: [1 x 3]: coefficient of planar fit
%
%
% [1] Wilczak, J. M., Oncley, S. P., & Stage, S. A. (2001).
% Sonic anemometer tilt correction algorithms.
% Boundary-Layer Meteorology, 99(1), 127-150.
%
% Author. E. Cheynet - University of Stavanger
%
%  last modified: 20.01.2018

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('method','PF');
p.addOptional('Err',1);
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
method = p.Results.method ;
Err = p.Results.Err ;

%% WILCZAK method is PF

[M,N]=size(u);
if M==1 && N ==1,
    error('meanU must be a vector, not  a scalar');
end


% If w is an empty variable, only the rotation in the horizontal plane is
% used. The additional parameter "method" is useless in this case
if isempty(w) %
    U1 = nan(size(u));
    V1 = nan(size(u));
    for ii=1:M
        
        A0 = [u(ii,:);v(ii,:)];
        R = nanmean(atan2(A0(1,:),A0(2,:)));
        R1 = [cos(R),-sin(R);sin(R),cos(R)];
        A1 = R1*A0;
        
        U1(ii,:) = A1(2,:);
        V1(ii,:) = A1(1,:);
    end
    W1= [];
    return
end

if strcmpi(method,'PF')
    meanU = nanmean(u,2)';
    meanV = nanmean(v,2)';
    meanW = nanmean(w,2)';
    meanU(isnan(meanU)) = [];
    meanW(isnan(meanW)) = [];
    meanV(isnan(meanV)) = [];
    M = numel(meanU); % updat the value of M
    U1 = nan(size(u));
    V1 = nan(size(u));
    W1 = nan(size(u));
    
    [b0,b1,b2,r1] = findB(meanU,meanV,meanW,M);
    
    
    if Err==0
        b0 = 0;
    elseif Err~=1
        error(' ''Err'' should be  0 or 1 ');
    end
    
    
    Deno = sqrt(1+b1.^2+b2.^2);
    p31 =-b1./Deno;
    p32 =-b2./Deno;
    p33 =1.00/Deno;
    
    cosGamma = p33./sqrt(p32.^2+p33.^2);
    sinGamma = -p32./sqrt(p32.^2+p33.^2);
    cosBeta = sqrt(p32.^2+p33.^2);
    sinBeta = p31;
    
    
    R2 = [1,0,0; 0,cosGamma,-sinGamma;0,sinGamma,cosGamma];
    R3 = [cosBeta,0,sinBeta;0,1,0;-sinBeta,0,cosBeta];
    
    A0 = R3'*R2'*[meanU;meanV;meanW];
    Alpha = atan2(A0(2,:),A0(1,:));
    
    
    for ii=1:M
        R1 = [cos(Alpha(ii)),-sin(Alpha(ii)),0;sin(Alpha(ii)),cos(Alpha(ii)),0;0,0,1];
        A1 = R1'*((R3'*R2')*[u(ii,:);v(ii,:);w(ii,:)-b0]);
        U1(ii,:) = A1(1,:);
        V1(ii,:) = A1(2,:);
        W1(ii,:) = A1(3,:);
        
    end
    
    if nargout ==4
        varargout{1} = [b0,b1,b2,r1];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif  strcmpi(method,'rot2')
    U1 = nan(size(u));
    V1 = nan(size(u));
    W1 = nan(size(u));
    for ii=1:M
        
        % first rotation (around z axis)
        A01 = [u(ii,:);v(ii,:)];
        R1 = atan2(nanmean(A01(2,:)),nanmean(A01(1,:)));
        %R1 = [cos(R1),-sin(R1);sin(R1),cos(R1)];
        R1 = [cos(R1),sin(R1);-sin(R1),cos(R1)];   %inversion de signe
        %dans les sinus?
        A1 = R1*A01;
        u1 = A1(1,:);
        V1(ii,:) = A1(2,:);
        
        % second rotation (around y axis)
        A02 = [u1;w(ii,:)];
        R2 = atan2(nanmean(A02(2,:)),nanmean(A02(1,:)));
        RotY = [cos(R2),sin(R2);-sin(R2),cos(R2)];
        
        A2 = RotY*A02;
        
        U1(ii,:) = A2(1,:);
        W1(ii,:) = A2(2,:);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif  strcmpi(method,'rot3')
    
    U1 = nan(size(u));
    V1 = nan(size(u));
    W1 = nan(size(u));
    for ii=1:M
        
        % first rotation (around z axis)
        A01 = [u(ii,:);v(ii,:)];
        R1 = atan2(nanmean(A01(2,:)),nanmean(A01(1,:)));
        R1 = [cos(R1),sin(R1);-sin(R1),cos(R1)];
        A1 = R1*A01;
        u1 = A1(1,:);
        v1 = A1(2,:);
        
        % second rotation (around y axis)
        A02 = [u1;w(ii,:)];
        R2 = atan2(nanmean(A02(2,:)),nanmean(A02(1,:)));
        RotY = [cos(R2),sin(R2);-sin(R2),cos(R2)];
        
        A2 = RotY*A02;
        
        U1(ii,:) = A2(1,:);
        w1 = A2(2,:);
        
        
        % third rotation (around x axis)
        A03 = [v1;w1];
        covVW = nanmean(v1(:).*w1(:));
        diffVW = nanvar(v1)-nanvar(w1);
        R3  = 0.5*atan2(2*covVW,diffVW);
        RotX = [cos(R3),sin(R3);-sin(R3),cos(R3)];
        
        A3 = RotX*A03;
        V1(ii,:) = A3(1,:);
        W1(ii,:) = A3(2,:);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error(' ''method'' must be ''PF'', ''rot2'' or ''rot3'' ')
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [b0,b1,b2,r1] = findB(meanU,meanV,meanW,M)
        % Code taken from [1] after correciton of the typos
        su=nansum(meanU);
        sv=nansum(meanV);
        sw=nansum(meanW);
        suv=meanU*meanV';
        suw=meanU*meanW';
        svw=meanV*meanW';
        su2=meanU*meanU';
        sv2=meanV*meanV';
        H=[M su sv; su su2 suv; sv suv sv2];
        g=[sw suw svw]'  ;
        x=H\g    ;
        r1=g-H*x;
        fprintf('%10.5f %10.5f %10.5f\n %6.5f %6.5f %6.5f',x,r1);
        fprintf('\n');
        
        
        b0 = x(1);
        b1 = x(2);
        b2 = x(3);
    end
end


