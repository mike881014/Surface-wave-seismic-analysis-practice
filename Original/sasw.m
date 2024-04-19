function [fe, Vphe]=sasw(r1, r2, dt, D, Lmin, Lmax)
% This function peforms the SASW analysis.
%
% [fe,Vphe]= SASW(r1,r2,dt,D,Lmin,Lmax)
%		fe = frequency
%		Vphe = corresponding phase velocity
%		r1 = seismic data of first receiver
%		r2 = seismic data of second receiver
%		dt = sampling rate of seismic data
% 		D = Geophone spacing
%		Lmin = Lmin*D is the minimum effective wavelength
%		Lmax = Lmax*D is the maximum effective wavelength

% Note that phase unwrapping in SASW is tricky

% Written by C-P Lin, Last modified 8/5/2002.

% Corss-Power Spectra
% R1(f) and R2(f) correspond to the Fourier transforms of time records from two receivers located a distance D apart.
% The bar above r1,r2(f) corresponds to the frequancy-domain average of several records.
% Parameter n is the number of records averaged.
% The asterisk above R2(f) corresponds to the complex conjugate operator.

% Time/Frequency
N=length(r1);
df=1/(N*dt);
f=[0:df:(N/2)*df]';

% FFT of Signals
R1f=fft(r1,N); 
R1f_conjugate=conj((R1f));
R2f=fft(r2,N); 
R2f_conjugate=conj((R2f));

psd1=(R1f.*R1f_conjugate); % Power Spectral Density (Autopower spectrum)
psd2=(R2f.*R2f_conjugate);
csd=(R1f.*R2f_conjugate); % Cross Spectral Density (Cross Power Spectrum)
%csd=(R2f.*R1f_conjugate); 

% Coherence Function
gamma_square=csd.*conj(csd)./(psd1.*psd2);

% phase shift (phi) unwraped
Nf=round(200/df); 
f=f(2:Nf+1); 
phi=angle(csd(2:Nf+1));
punwrap = 1.0; % parameter for unwrapping
phi_un=unwrap(phi,punwrap*pi)*180/pi; 

% travel time (t)
TravelTime=phi_un./(360*f);

% phase velocity (Vph)
% distance between the receivers (D)
Vph=D./TravelTime;

% wavelength (Lph)
Lph=Vph./f;

% Frequency range (filtering criteria)
fe=[]; Vphe=[];
for i=1:length(f)
   if Lph(i)>(Lmin*D) & Lph(i)<(Lmax*D)
      fe=[fe; f(i)];
      Vphe=[Vphe; Vph(i)];
   end
end

% Plotting of phase difference as a function of frequency
figure
plot(f,phi,'b','LineWidth',2); hold on;
set(gca, 'LineWidth',2, 'FontSize',16');
xlabel('f(Hz)');
ylabel('\phi (rad)');
set(gca, 'LineWidth',2, 'FontSize',16');
plot(f,phi_un*pi/180,'r','LineWidth', 2);
xlabel('frequency(Hz)');
ylabel('\Delta\phi (rad)');