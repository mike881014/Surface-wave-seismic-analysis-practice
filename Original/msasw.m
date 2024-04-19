function [f,Vph, R2, UangleU, Uangle]=MSASW(Data,dt,x,fmin,fmax,df)
% This function peforms the multi-channel SASW analysis.
%
% [f,Vph]=MSASW(Data,dt,dx,vmin,vmax,Nv,fmin,fmax)
%		f = frequency
%		Vph = Corresponding phase velocity
%		Data = seismic data (No. of samples, No. of Channel)
%		dt = sampling rate of seismic data
% 		x = geophone position in column form
% 		fmin, fmax = frequency range
%		df = frequency resolution
%
% Written by C-P Lin, Last modified 8/5/2002.

[N, Ch]=size(Data); %No. of data, No. of channel

% Set up frequency resolution
Nf=ceil(1/(df*dt));
df=1/(Nf*dt);
Nf1=floor(fmin/df)+1;
Nf2=ceil(fmax/df)+1;
f=[(Nf1-1)*df:df:(Nf2-1)*df]'; % Frequency Range
%x=[0:dx:(Ch-1)*dx]'; % Offsets

U=fft(Data,Nf); % Fast Fourier Transformation
U=U(Nf1:Nf2,:);
Uangle=angle(U);
p_unwrap = 1.0;
UangleU=unwrap(angle(U),[p_unwrap*pi],2); % Unwrap phase angle



Vph=zeros(length(f),1);
R2=zeros(length(f),1);
for i=1:length(f)
   phi=UangleU(i,:);
   a=1; b=Ch; % Effective channel from a to b
   %b=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis
   [b bint r rint stats]=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis with R2 statistic
   Vph(i)=abs(b(1))*2*pi*f(i);
   R2(i)=stats(1); % R2 statistic
end

% Plotting of dispersion curve and R^2 as a function of frequency
%figure(5)
%subplot(211); hold on;
%plot(f,Vph,'b:');
%figure(5)
%subplot(212)
%plot(f,R2,'r:'); hold on;

% Plotting of phi(x) for a particular frequency
nf=36;
figure;
plot(x,Uangle(nf,:),'ro-');hold on;
plot(x,UangleU(nf,:),'LineWidth',2);hold on;
set(gca,'FontSize',16,'LineWidth',2);
xlabel('Offset (m)');
ylabel('\phi (rad)');
