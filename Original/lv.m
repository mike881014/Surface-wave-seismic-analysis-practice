function [v,l,A,U]=lf(Data,dt,x,vmin,vmax,dv,lmin,lmax,dl)
% This function peforms the lamda-frequency analysis of multi-channel surface wave by Lin(2004)
% This function includes windowing in time domain and space domain.  
% It also accounts for relative amplitude of different frequency in a correct manner.
%
% [f,l,A,U]=lf(Data,dt,x,fmin,fmax,df,lmin,lmax,dl)
%		f = frequency range
%		l = wavelength range
%		A = amplitude spectrum of overtone analysis
%		Data = seismic data (No. of samples, No. of Channel)
%		dt = sampling rate of seismic data
% 		x = geophone position in column form
% 		fmin, fmax = frequency range
%		df = frequency resolution	
% 		lmin, lmax specify wavelength range
%      dl = wavelength resolution
%
% Written by C-P Lin, Last modified 6/14/2004.

[N, Ch]=size(Data); %No. of data, No. of channel

% Trace Scaling (Energy Balancing)
for i=1:Ch
   sf=sqrt(sum(Data(:,i).^2)/N);
   Data(:,i)=Data(:,i)/sf;
end

% apply time window when vibroseis source is used
%wt=hamming(N);
%for i=1:Ch
%   Data(:,i)=Data(:,i).*wt;
%end

% Set up frequency and velocity resolution
t=0:dt:(N-1)*dt;
Nl=ceil((lmax-lmin)/dl);
l=[lmin:dl:lmin+(Nl-1)*dl]'; % Wavelength Range
Nv=ceil((vmax-vmin)/dv);
v=[vmin:dv:vmin+(Nv-1)*dv]'; % Velocity Range

% Windowing in space domain
W=[];
%w=bartlett(Ch)';
%w=hamming(Ch)';
w=kaiser(Ch,4)';
%w=ones(1,Ch);
for i=1:length(t)
   W=[W; w];
end

% Discrete space Fourier transform

% Discrete space Fourier transform
j=sqrt(-1);
for i=1:length(l) % Phase Shift (phi)
   [xx, tt]=meshgrid(x,t);
   tmp=exp(-j*2*pi/l(i)*xx); % phi=w/v=2*pi*f/v
   V=(Data.*W).*tmp; % Integral Transformation (with windowing)
   U(:,i)=sum(V')'; % Summing up
end


% Fast Fourier Transformation in Time Domain
for i=1:length(v)
    [ll,tt]=meshgrid(l,t);
    tmp=exp(j*2*pi*v(i)./ll.*tt);
    UV=U.*tmp;
    UU(i,:)=sum(UV);
end

A=abs(UU); % Amplidtude Spectrum 