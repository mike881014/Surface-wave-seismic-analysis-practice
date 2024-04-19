clear all,close all,clc
cd('/Users/gaoqianhui/Desktop/模擬/學姊');
% In the practicing now:
% Z: vertical direction, 
% X: radial direction, 
% Y: transversal direction

%% Simulated files combine
first_fnu = 1; last_fnu = 20; % The start and end station for analysis
NF=last_fnu-first_fnu+1; % Number of files

fn1='DB.X';
Direction_X='X';fn3_X=['.FX',Direction_X,'.semd'];
Direction_Y='Y';fn3_Y=['.FX',Direction_Y,'.semd'];
Direction_Z='Z';fn3_Z=['.FX',Direction_Z,'.semd'];

for i=1:1:NF
    fn2=num2str(first_fnu+i-1); % The following lines putting the pieces of file names together
    FileName_X=[fn1 fn2 fn3_X];
    FileName_Y=[fn1 fn2 fn3_Y];
    FileName_Z=[fn1 fn2 fn3_Z];

    RawData_X=load(FileName_X);
    RawData_Y=load(FileName_Y);
    RawData_Z=load(FileName_Z);

    y_X(:,i)=-RawData_X(:,2); % Add negitive sign for direction correction
    y_X_norm(:,i) = y_X(:,i)/max(abs(y_X(:,i)));
    y_Y(:,i)=-RawData_Y(:,2); % Add negitive sign for direction correction
    y_Y_norm(:,i) = y_Y(:,i)/max(abs(y_Y(:,i)));
    y_Z(:,i)=RawData_Z(:,2);
    y_Z_norm(:,i) = y_Z(:,i)/max(abs(y_Z(:,i)));

end

[N, Ch]=size(y_Z);
% The offset for analysis (Carefull)
dx=0.1; % Geophone spacing
x0=0.5; % distance between source & 1st receiver
x=[x0:dx:x0+(Ch-1)*dx]';

t = RawData_X(:,1); dt=t(2)-t(1); % time step
% [m,t0_ind] = min(abs(t)); % Time correction, the simulation time start from negitive time 
% t = t(t0_ind:end); N = length(t);
% y_X = y_X(t0_ind:end,:); y_Y = y_Y(t0_ind:end,:); y_Z = y_Z(t0_ind:end,:);
N = length(t); t=[0:dt:(N-1)*dt];

%% Figure plot in time domain
figure,
subplot(131),wigb(y_Z_norm,1,x,t);
set(gca,'FontSize',10,'Ydir','reverse');ylabel('Time (s)'); xlabel('Offset (m)'); title('Vertical component')
subplot(132),wigb(y_X_norm,1,x,t);
set(gca,'FontSize',10,'Ydir','reverse');ylabel('Time (s)'); xlabel('Offset (m)'); title('Radial component')
subplot(133),wigb(y_Y_norm,1,x,t);
set(gca,'FontSize',10,'Ydir','reverse');ylabel('Time (s)'); xlabel('Offset (m)'); title('Transversal component')

%% Zero padding to 1 sec
N = round(1/dt);t = [0:dt:(N-1)*dt];
data_length = length(y_X);
y_X = [y_X;zeros(N-data_length,length(x))];
y_Y = [y_Y;zeros(N-data_length,length(x))];
y_Z = [y_Z;zeros(N-data_length,length(x))];


%% Data muting
% y_mute = seismic_data_muting(x,y_Z,dt,N);
% % plot muted data
% nData=y_mute; 
% for i = 1:length(x)
%     nData(:,i) = y_mute(:,i)./max(abs(y_mute(:,i)));
% end
% figure;
% wigb(nData,1,x,t);
% ylabel('Time, sec'); xlabel('Source-receiver Distance, m'); set(gca,'fontsize',14)

%% Dispersion analysis
% Analysis frequency range
fmin=0; fmax=14000; df=2000;
% Analysis velocity range
vmin = 0; vmax=2000; dv=500;
% Analysis wave length range
lmin=1; lmax=100; dl=0.5;

% Peform f-v transform, f-v

[f,v,A_Z,Z_fx]=fv(y_Z,dt,x,fmin,fmax,df,vmin,vmax,dv);
[f,v,A_X,X_fx]=fv(y_X,dt,x,fmin,fmax,df,vmin,vmax,dv);
[f,v,A_Y,Y_fx]=fv(y_Y,dt,x,fmin,fmax,df,vmin,vmax,dv);

% Perform f-l transform and then transform to l-v
[f,l,Al_Z]=fl(y_Z,dt,x,fmin,fmax,df,lmin,lmax,dl);[l,v, A_lv_Z]=fl_to_lv(Al_Z, f, l, v);
[f,l,Al_X]=fl(y_X,dt,x,fmin,fmax,df,lmin,lmax,dl);[l,v, A_lv_X]=fl_to_lv(Al_X, f, l, v);
[f,l,Al_Y]=fl(y_Y,dt,x,fmin,fmax,df,lmin,lmax,dl);[l,v, A_lv_Y]=fl_to_lv(Al_Y, f, l, v);

fe=f;
Vphe_Z=zeros(length(fe),1);Vphe_X=zeros(length(fe),1);Vphe_Y=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(A_Z(i,:));
   Vphe_Z(i)=v(Ai);
   An_Z(i,:) = A_Z(i,:)/Amax;

   [Amax Ai]=max(A_X(i,:));
   Vphe_X(i)=v(Ai);
   An_X(i,:) = A_X(i,:)/Amax;

   [Amax Ai]=max(A_Y(i,:));
   Vphe_Y(i)=v(Ai);
   An_Y(i,:) = A_Y(i,:)/Amax;
end


%% Radial to vertical relation (sRVSR, RVPD and RVSR)
CWU=conj(Z_fx).*X_fx;
RVSR = abs(X_fx)./abs(Z_fx); 
RVPD=angle(CWU); 
sRVSR=sign(imag(CWU)).*RVSR;  

% Subplot: sRVSR, RVPD, RVSR
figure;
subplot(131),pcolor(x,f,sRVSR); 
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse');
colormap(gca,jet);colorbar;caxis([-2 2]);
title('sRVSR'); xlabel('x (m)'); ylabel('f (Hz)');

subplot(132),pcolor(x,f,RVPD);
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse');
colormap(gca,jet);colorbar;caxis([-pi pi]);
title('RVPD'); xlabel('x (m)'); ylabel('f (Hz)');

subplot(133),pcolor(x,f,RVSR);
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
colormap(gca,gray);colorbar;caxis([-2 2]);
title('RVSR'); xlabel('x (m)'); ylabel('f (Hz)');

%% Plot in f-x domain
figure,
subplot(131),pcolor(x,f,abs(Z_fx));
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
% colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 
colormap(gca,jet);colorbar;%caxis([-2 2]);
title('Vertical'); xlabel('x (m)'); ylabel('f (Hz)');

subplot(132),pcolor(x,f,abs(X_fx));
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
% colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 
colormap(gca,jet);colorbar;%caxis([-2 2]);
title('Radial'); xlabel('x (m)'); ylabel('f (Hz)');
 
subplot(133),pcolor(x,f,abs(Y_fx));
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
% colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 
colormap(gca,jet);colorbar;%caxis([-2 2]);
title('Transversal'); xlabel('x (m)'); ylabel('f (Hz)');


%% Plot in fv domain
figure;
subplot(131),imagesc(f,v,An_Z);colormap(jet);
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);
hold on;plot(Vphe_Z,fe,'w.','Markersize',40);
ylabel('Phase Velocity, V_{ph} (m/s)');xlabel('Frequency, f (Hz)'); title('Vertical')
set(gca,'fontsize',14);

% subplot(132),imagesc(v,f,An_X);colormap(jet);
% % cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);
% hold on;plot(Vphe_X,fe,'w.','Markersize',10);
% xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Frequency, f (Hz)'); title('Radial')
% set(gca,'fontsize',14);
% 
% subplot(133),imagesc(v,f,An_Y);colormap(jet);
% % cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);
% hold on;plot(Vphe_Y,fe,'w.','Markersize',10);
% xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Frequency, f (Hz)'); title('Transversal')
% set(gca,'fontsize',14);

%% Plot in l-v domain
figure;
subplot(131),imagesc(v,l,A_lv_Z);colormap(jet);
cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);caxis([0 1]);
hold on;plot(Vphe_Z,Vphe_Z./fe,'w.','Markersize',10);
xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Wavelength, \lambda (m)'); title('Vertical')
set(gca,'fontsize',14);

subplot(132),imagesc(v,l,A_lv_X);colormap(jet);
cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);caxis([0 1]);
hold on;plot(Vphe_X,Vphe_X./fe,'w.','Markersize',10);
xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Wavelength, \lambda (m)'); title('Radial')
set(gca,'fontsize',14);

subplot(133),imagesc(v,l,A_lv_Y);colormap(jet);
cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);caxis([0 1]);
hold on;plot(Vphe_Y,Vphe_Y./fe,'w.','Markersize',10);
xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Wavelength, \lambda (m)'); title('Transversal')
set(gca,'fontsize',14);


%% Convolution
% Source type: Half sine
% Tp=0.01;
% tp=[0:dt:Tp]'; Np=length(tp);
% Pt=sin(2*pi/(Tp*2)*tp)*1; % Half sine source
% N = length(t);
% Pt = [Pt;zeros(N-Np,1)];
% 
% y_Offset = y(:,find(x==Offset)); % Examine the time domain data at certain offset
% Single = conv(y_Offset,Pt);
% % Single = -Single/max(Single);
% figure,plot(t,-Single(1:length(t)),'r-'),xlabel('Time (s)'),title(['y*Pt at offset=',num2str(Offset)])
% set(gca,'fontsize',14);

%% try to mix f-k
% This is an example of Multichannel Analysis of Surface Wave
% using f-k analysis

[N, Ch]=size(Data);x1=10; x=[x1:dx:x1+(Ch-1)*dx]';
[f,k,A]=fk(Data,dt,x,fmin,fmax,df,kmin,kmax,dk);


% Find the dispersion curve from the f-k image
fe=f;
Vphe=zeros(length(fe),1);
ke=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(A(i,:));
   ke(i)=k(Ai);
   Vphe(i)=2*pi*fe(i)/k(Ai);
end

% Plotting the f-k image and dispersion curve
figure;
imagesc(k,f,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=1;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
xlabel('Wavenumber');
ylabel('Frequency, f (Hz)');
hold on;
plot(ke,fe,'r.');hold on;


% Plotting the f-v image and dispersion curve
[kk,ff]=meshgrid(k,f);
vv=(2*pi*ff)./kk;

figure
surfc(vv,ff,A);
view([0,-90]);
shading flat;

colormap(meshgrid(1:-1/63:0,1:3)');
%c1=0; c2=0; c3=0; cs=3;
%cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
%cmap=cmap.^cs;
%colormap(cmap);
colorbar;
%set(gca,'xticklabel',{v(20),v(40),v(60),v(80),v(100),v(120)});
%set(gca,'yticklabel',{f(10),f(20),f(30),f(40),f(50),f(60),f(70),f(80),f(90)});
title('Overtone Diagram');
xlabel('Phase Velocity, v_p_h (m/s)');
ylabel('Frequency, f (Hz)');
hold on;
plot(Vphe,fe,'r.');hold on;
axis([0 1000 0 80]);


% Perform f-v Analysis
vmin=150; vmax=500; dv=1; % Velocity Range
fmin=5; fmax=80; df=1; % Frequency Range
dx=1; % Geophone spacing
[N, Ch]=size(Data);x1=28; x=[x1:dx:x1+(Ch-1)*dx]';
[f,v,A1]=fv(Data,dt,x,fmin,fmax,df,vmin,vmax,dv);

% Find the dispersion curve from the f-v image
fe1=f;
Vphe1=zeros(length(fe1),1);
for i=1:length(fe1);
   [Amax Ai]=max(A1(i,:));
   Vphe1(i)=v(Ai);
end

% Perform mutichanel SASW analysis
%[fe2, Vphe2]=MSASW(Data,dt,x,fmin,fmax,df);

% Compare and plot f-k, f-v, and MSASW
figure
plot(fe, Vphe,'r*');hold on;
plot(fe1, Vphe1,'bo');

%plot(fe2, Vphe2,'b.');
%axis([0 80 0 500]);
ylabel('Phase Velocity, v_p_h (m/s)');
xlabel('Frequency, f (Hz)');
%legend('f-k','f-v','MSASW');


