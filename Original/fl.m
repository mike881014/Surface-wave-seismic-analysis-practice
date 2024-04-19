function [f,l,A,U]=fl(Data,dt,x,fmin,fmax,df,lmin,lmax,dl)
% This function peforms the overtone analysis of multi-channel surface wave by Lin(2002)
% This function includes windowing in time domain and space domain.  
% It also accounts for relative amplitude of different frequency in a correct manner.
%
% [f,v,A]=fl(Data,dt,x,fmin,fmax,df,lmin,lmax,dl)
%		f = frequency range
%		l = wavelength range
%		A = amplitude spectrum of overtone analysis
%		Data = seismic data (No. of samples, No. of Channel)
%		dt = sampling rate of seismic data
% 		x = geophone position in column form
% 		fmin, fmax = frequency range
%		df = frequency resolution	
% 		lmin, lmax specify wavelength range
%     dl = wavelength resolution
%
% Written by C-P Lin, Last modified 10/14/2003.

[N, Ch]=size(Data); %No. of data, No. of channel

Nf=ceil(1/(df*dt));
df=1/(Nf*dt);
Nf1=floor(fmin/df)+1;
Nf2=ceil(fmax/df)+1;
f=[(Nf1-1)*df:df:(Nf2-1)*df]'; % Frequency Range
Nl=ceil((lmax-lmin)/dl);
l=[lmin:dl:lmin+(Nl-1)*dl]'; % Velocity Range

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

% Fast Fourier Transformation in Time Domain
U=fft(Data,Nf); % Fast Fourier Transformation
U=U(Nf1:Nf2,:);

% Windowing in space domain
W=[];
%w=bartlett(Ch)';
%w=hamming(Ch)';
w=kaiser(Ch,4)';
%w=ones(1,Ch);
for i=1:length(f)
   W=[W; w];
end

% Discrete space Fourier transform
j=sqrt(-1);
for i=1:length(l) % Phase Shift (phi)
   [xx, ff]=meshgrid(x,f);
   tmp=exp(j*2*pi/l(i).*xx); % phi=w/v=2*pi*f/v
   V=(U.*W).*tmp; % Integral Transformation (with windowing)
   UU(:,i)=sum(V')'; % Summing up
end
A=abs(UU); % Amplidtude Spectrum
