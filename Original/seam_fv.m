% f-v results for pseudo-section survey line with seaming
function [f,v,A,Y,y,YU] = seam_fv(y,dt,xY,fmin,fmax,df,vmin,vmax,dv)

% This function seaming the phase difference at overlapped offsets
%
% Input:
%       y = data, m x n x p
%       dt = sampling rate, sec
%       xY = space locatoin for each geophone after seaming, m
%       fmin = min frequency, Hz
%       fmax = max frequency, Hz
%       df = frequency space, Hz
%       vmin = min interested velocity, m/s
%       vmax = max interested velocity, m/s
%       dv = velocity space, m/s
% Output:
%       f = frequency series, Hz
%       v = phase velocity, m/s
%       A = spectrum amplitudein f-k domain
%       Y = Spectrum in frequency domain
% 		y = energy balance data y(t,x)
%       YU = unseaming spectrum in frequency domain, m x n x p
%
% Written by C.H, Lin, Last modified Mar. 2009.

[M N NF]=size(y);

%% Trace Scaling (Energy Balancing)
for i=1:NF
	for k=1:N
   	sf=sqrt(sum(y(:,k,i).^2)/M);
	   y(:,k,i)=y(:,k,i)/sf;
   end
end

%% Set up frequency and velocity resolution
Nf=ceil(1/(df*dt));
df=1/(Nf*dt);
Nf1=floor(fmin/df)+1;
Nf2=ceil(fmax/df)+1;
f=[(Nf1-1)*df:df:(Nf2-1)*df]'; % Frequency Range
Nv=ceil((vmax-vmin)/dv);
v=[vmin:dv:vmin+(Nv-1)*dv]'; % Velocity Range


%% Seaming in frequency domain
Y = fft(y,Nf);
Y = Y(Nf1:Nf2,:,:);
Yangle = angle(Y);
[M N NF] = size(Y);
YU = Y;

% correct Yangle, Y, and y
j=sqrt(-1);
for i=1:NF-1
   dangle=Yangle(:,N,i)-Yangle(:,1,i+1);
   dangle=repmat(dangle,1,N);
   shift=exp(j*dangle);
   Y(:,:,i+1)=Y(:,:,i+1).*shift;
   Yangle(:,:,i+1)=angle(Y(:,:,i+1));
end
Y=[reshape(Y(:,1:N-1,:),M,(N-1)*NF) Y(:,N,NF)];

[M Ch]=size(Y);


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
for i=1:length(v) % Phase Shift (phi)
   [xx, ff]=meshgrid(xY,f);
   tmp=exp(j*2*pi*ff/v(i).*xx); % phi=w/v=2*pi*f/v
   V=(Y.*W).*tmp; % Integral Transformation (with windowing)
   UU(:,i)=sum(V')'; % Summing up
end
A=abs(UU); % Amplitude Spectrum


