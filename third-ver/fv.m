function [f,v,A,U]=fv(Data,dt,x,fmin,fmax,df,vmin,vmax,dv)
% This function peforms the f-v analysis of multi-channel surface wave by Lin(2004)
% This function includes windowing in time domain and space domain.  
% It also accounts for relative amplitude of different frequency in a correct manner.
%
% [f,v,A]=fv(Data,dt,x,fmin,fmax,df,vmin,vmax,dv)
%		f = frequency range
%		v = velocity range
%		A = amplitude spectrum of overtone analysis
%		Data = seismic data (No. of samples, No. of Channel)
%		dt = sampling rate of seismic data
% 		x = geophone position in column form
% 		fmin, fmax = frequency range
%		df = frequency resolution	
% 		vmin, vmax specify velocity range
%       dv = velocity resolution
%
% Written by C-P Lin, Last modified 6/19/2005.

[N, Ch]=size(Data); %No. of data, No. of channel

% Trace Scaling (Energy Balancing)
% for i=1:Ch
%    sf=sqrt(sum(Data(:,i).^2)/N);
%    Data(:,i)=Data(:,i)/sf;
% end

% Apply time window when necessary (e.g. vibroseis source)
%Wt=[];
%wt=bartlett(N);
%wt=hamming(N);
%wt=kaiser(N,4);
%wt=ones(N,1);
%for i=1:Ch
%   Wt=[Wt wt];
%end
%Data=Data.*Wt; 

% Set up frequency and velocity resolution
Nf=ceil(1/(df*dt));
df=1/(Nf*dt);
Nf1=floor(fmin/df)+1;
Nf2=ceil(fmax/df)+1;
f=[(Nf1-1)*df:df:(Nf2-1)*df]'; % Frequency Range
Nv=ceil((vmax-vmin)/dv);
v=[vmin:dv:vmin+(Nv-1)*dv]'; % Velocity Range

% Fast Fourier Transformation in Time Domain
U=fft(Data,Nf); % Fast Fourier Transformation
U=U(Nf1:Nf2,:);
% Windowing in space domain
W=[];
%w=bartlett(Ch)';
%w=hamming(Ch)';
% w=kaiser(Ch,4)';
w=ones(1,Ch);
for i=1:length(f)
   W=[W; w];
end
Uf=U.*W;

% Discrete space Fourier transform --> velocity domain
j=sqrt(-1);
for i=1:length(v) 
   [xx, ff]=meshgrid(x,f);
   tmp=exp(j*2*pi*ff/v(i).*xx); 
   V=(Uf).*tmp; 
   UU(:,i)=sum(V')'; 
end
A=abs(UU); % f-v Spectrum
