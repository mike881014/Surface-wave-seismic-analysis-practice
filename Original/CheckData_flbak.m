% This is an example of Multichannel Analysis of Surface Wave
% using overtone analysis
% Read seismic data
[y, dt]=read_org('Mc403_SX.org');
dx=1; % Geophone spacing

%y=fliplr(y);
%y=y(:,4:24);

%---------------------------------------------------------------------------------------
% Peform Overtone Analysis
lmin=1; lmax=40; dl=0.1; % Wavelength Range
fmin=1; fmax=100; df=1; % Frequency Range

[N, Ch]=size(y);
x1=15; x=[x1:dx:x1+(Ch-1)*dx]';
t=0:dt:(N-1)*dt;

figure
wigb(y,1.,x,t);
set(gca,'FontSize',10);
ylabel('Time (s)')
xlabel('Offset (m)')
set(gca,'XAxisLocation','top','Ydir','reverse');

[f,l,A]=fl(y,dt,x,fmin,fmax,df,lmin,lmax,dl);

% Find the dispersion curve from the f-l image
fe=f;
Vphe=zeros(length(fe),1);
le=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(A(i,:));
   le(i)=l(Ai);
   Vphe(i)=fe(i)*l(Ai);
end

%nf=50;
%[maxA, IA]=max(A(nf,:));
%plot(v,A(nf,:)/maxA*f(nf));
%plot(v(IA),f(nf),'b^');
% ------------------------------------------------------------------------------------
Nf=ceil(1/(df*dt));
df=1/(Nf*dt);
Nf1=floor(fmin/df)+1;
Nf2=ceil(fmax/df)+1;
f=[(Nf1-1)*df:df:(Nf2-1)*df]'; % Frequency Range

U=fft(y,Nf); % Fast Fourier Transformation
U=U(Nf1:Nf2,:);
p_unwrap = 1.0;
PhiU=unwrap(angle(U),[p_unwrap*pi],2); % Unwrap phase angle

%figure
%for i=5:10:95;
%   tmp=10*(PhiU(i,:)-PhiU(i,1))/(PhiU(i,1)-PhiU(i,Ch))+f(i);
%   plot(x, tmp); hold on;
%end



% Perform mutichanel SASW analysis
[fe2, Vphe2, R2, PhiU, Phi]=MSASW(y,dt,x,fmin,fmax,df);
%plot(Vphe2,fe2,'bo');
%legend('MWTSW','MSASW');

figure
subplot(131);
wigb(abs(U),1.,x,f);
set(gca,'FontSize',10);
ylabel('frequency (Hz)')
xlabel('Offset (m)')
set(gca,'Ydir','reverse');


subplot(132);
plot(R2,fe2,'b'); hold on;
axis([0 1.0 0 100]);
set(gca,'Ydir','reverse','FontSize',10);
ylabel('Frequency (Hz)','FontSize',10);
xlabel('R^2','FontSize',10);

% Plotting the overtone image and dispersion curve
subplot(133);
imagesc(l,f,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=3;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
%set(gca,'xticklabel',{v(20),v(40),v(60),v(80),v(100),v(120)});
%set(gca,'yticklabel',{f(10),f(20),f(30),f(40),f(50),f(60),f(70),f(80),f(90)});
%title('Overtone Diagram');
set(gca,'FontSize',10);
xlabel('Wavelength (m)');
ylabel('Frequency (Hz)');
hold on;
plot(le(4:length(le)-2),fe(4:length(le)-2),'r.');hold on;