clear all;clc%close all;clc
cd("/Users/gaoqianhui/Desktop/模擬/學姊")
% cd('C:\Users\Tsai Jung\Dropbox\SPEC3D\vs=200_vp=1428_xyz=2008050');
% b=0.5;ht=2.5;h=0.5;
% hold on;
% plot(b,-ht*ones(length(b),1)-h,'-');
% figure,plot(b,-ht*ones(length(b),1),'-k');
% hold on;
% plot(b,-ht*ones(length(b),1)-h,'-k');
% plot([10:1:50],zeros(length([10:1:50]),1),'r^');
% set(gcf,'color','w');
% figure;
% subplot(141)
% line([29.75,30.25],[-ht,-ht]);
% subplot(142)
% line([29.75,30.25],[-ht-h,-ht-h]);
% subplot(143)
% line([29.75,29.75],[-ht,-ht-h]);
% subplot(144)
% line([30.25,30.25],[-ht,-ht-h]);
% angles = linspace(0, 2*pi, 500);
% radius = 1;
% CenterX = 30;
% CenterY = -5;
% x = radius * cos(angles) + CenterX;
% y = radius * sin(angles) + CenterY;
% plot(x, y, 'b');
% hold on;
% plot(CenterX, CenterY, 'k+');
% grid on;
% axis equal;
% xlabel('X', 'FontSize', 20);
% ylabel('Y', 'FontSize', 20);
first_fnu = 1;
last_fnu =  20;
NF=last_fnu-first_fnu+1; % Number of files
fn1='DB.X';Direction1='Z';Offset = 50;
fn3_Z=['.FX',Direction1,'.semd'];

Direction2='X';
fn3_X=['.FX',Direction2,'.semd'];

Direction3='Y';
fn3_Y=['.FX',Direction3,'.semd'];

for i=1:1:NF
    fn2=num2str(first_fnu+i-1); % The following lines putting the pieces of file names together
    FileName1=[fn1 fn2 fn3_Z];
    RawData=load(FileName1);
    yZ(:,i)=RawData(:,2);
    
    FileName2=[fn1 fn2 fn3_X];
    RawData=load(FileName2);
    yX(:,i)=-RawData(:,2);
    
    FileName3=[fn1 fn2 fn3_Y];
    RawData=load(FileName3);
    yY(:,i)=-RawData(:,2);

    %y(:,i) = smooth(y(:,i),10);
    yZ_norm(:,i) = yZ(:,i)/max(abs(yZ(:,i)));
    yX_norm(:,i) = yX(:,i)/max(abs(yX(:,i)));
    yY_norm(:,i) = yY(:,i)/max(abs(yY(:,i)));
end



[N, Ch]=size(yZ);
dx=0.1; % Geophone spacing
x1=0.5; % distance between source & 1st receiver
y1=0.2;
x=[x1:dx:x1+(Ch-1)*dx]';
y=[y1:dx:y1+(Ch-1)*dx]';
% t=0:dt:(N-1)*dt;
t = RawData(:,1); dt=t(2)-t(1); % time step
% Dt=50;
[m,t0_ind] = min(abs(t)); % Time correction, the simulation time start from negitive time 
t = t(t0_ind:end);
yZ = yZ(t0_ind:end,:);
yX = yX(t0_ind:end,:);
yY = yY(t0_ind:end,:);
N = length(t);


%% Figure plot
figure;
subplot(121)
wigb(yZ,0.1,x,t);
set(gca,'FontSize',18);
set(gcf,'color','w');
ylabel('Time (s)')
xlabel('Offset (m)')
if Direction1=='X',
    title('Radial component')
elseif Direction1=='Z',
    title('Vertical component')
end
set(gca,'Ydir','reverse');
% legend('(x,y,z)=(200,80,25),mesh size=0.8m')
subplot(122)
wigb(yX,1.,x,t);
set(gca,'FontSize',18);
set(gcf,'color','w');
ylabel('Time (s)')
xlabel('Offset (m)')
if Direction2=='X',
    title('Radial component')
elseif Direction2=='Z',
    title('Vertical component')
end
set(gca,'Ydir','reverse');
% legend('(x,y,z)=(200,80,25),mesh size=0.8m')
N=round(1/dt);t = [0:dt:(N-1)*dt];
data_length = length(yX);
yX=[yX;zeros(N-data_length,length(x))];
yY=[yY;zeros(N-data_length,length(x))];
yZ=[yZ;zeros(N-data_length,length(x))];
% %% Plot time domain data at vertain offset
% y_Offset = yZ(:,find(x==Offset)); % Examine the time domain data at certain offset
% figure,plot(t,y_Offset,'r-'),xlabel('Time (s)'),title(['Time domain waveform at offset=',num2str(Offset)]);legend('(x,y,z)=(200,80,25),mesh size=0.8m')
% set(gcf,'color','w');
% set(gca,'fontsize',30);

% %% Data muting
% yZ = seismic_data_muting(x,yZ,dt,N);
% % plot muted data
% nData=yZ; 
% for i = 1:length(x)
%     nData(:,i) = yZ(:,i)./max(abs(yZ(:,i)));
% end
% figure;
% wigb(nData,1,x,t);
% ylabel('Time, sec'); xlabel('Source-receiver Distance, m'); set(gca,'fontsize',14)


% Analysis frequency range
fmin=5; fmax=80; df=0.5;
% Analysis velocity range
vmin=100; vmax=500; dv=0.5;
% Analysis wave length range
lmin=0.5; lmax=20; dl=0.1;

 %% Peform f-v transform, f-v
[f,v,Av_Y,U_fY]=fv(yY,dt,x,fmin,fmax,df,vmin,vmax,dv);
[f,v,Av_X,U_fX]=fv(yX,dt,x,fmin,fmax,df,vmin,vmax,dv);
[f,v,Av_Z,U_fZ]=fv(yZ,dt,x,fmin,fmax,df,vmin,vmax,dv);
fe=f;
Vphe=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(Av_Z(i,:));
   Vphe_Z(i)=v(Ai);
   Avn_Z(i,:) = Av_Z(i,:)/Amax;

   [Amax Ai]=max(Av_X(i,:));
   Vphe_X(i)=v(Ai);
   Avn_X(i,:) = Av_X(i,:)/Amax;

   [Amax Ai]=max(Av_Y(i,:));
   Vphe_Y(i)=v(Ai);
   Avn_Y(i,:) = Av_Y(i,:)/Amax;

end

 

% figure;
% 
% imagesc(v,f,Avn_Z);
% % colormap(meshgrid(1:-1/63:0,1:3)'); 
% colormap(jet)
% hold on;plot(Vphe_Z(1:4:152),fe(1:4:152),'w.','Markersize',20);
% shading interp;
% title('Vertical offset=20~30m')
% xlabel('Phase Velocity, v_{ph} (m/s)');ylabel('Frequency, f (Hz)');
% axis([100 500 5 80]);
% %% Plot in fv domain
% 
% figure;
% subplot(121)
% imagesc(v,f,Avn_Z);
% % colormap(meshgrid(1:-1/63:0,1:3)'); 
% colormap(jet)
% hold on;plot(Vphe_Z(1:4:152),fe(1:4:152),'w.','Markersize',20);
% shading interp;
% title('Vertical offset=20~30m')
% xlabel('Phase Velocity, v_{ph} (m/s)');ylabel('Frequency, f (Hz)');
% axis([100 500 5 80]);
% 
% % title('Dispersion curve')
% set(gcf,'color','w');
% set(gca,'fontsize',14);
% subplot(122)
% shading interp;
% imagesc(v,f,Avn_X);
% % colormap(meshgrid(1:-1/63:0,1:3)'); 
% colormap(jet)
% hold on;plot(Vphe_X(1:4:152),fe(1:4:152),'w.','Markersize',20);
% shading interp;
% title('Radial offset=')
% xlabel('Phase Velocity, v_{ph} (m/s)');ylabel('Frequency, f (Hz)');
% axis([100 500 5 80]);
% % title('Dispersion curve')
% set(gcf,'color','w');
% set(gca,'fontsize',14);

% Windowed Plot in fv domain

% figure;
% subplot(2,5,1)
% imagesc(v,f,Avn_Z);
% % colormap(meshgrid(1:-1/63:0,1:3)'); 
% colormap(jet)
% hold on;plot(Vphe_Z(1:4:182),fe(1:4:182),'w.','Markersize',20);
% shading interp;
% 
% title(['Vertical vc=',num2str(last_fnu)])
% xlabel('Phase Velocity, v_{ph} (m/s)');ylabel('Frequency, f (Hz)');
% axis([100 350 10 100]);
% % title('Dispersion curve')
% set(gcf,'color','w');
% set(gca,'fontsize',18);

% figure;
% subplot(2,5,2)
% shading interp;
% imagesc(v,f,Avn_X);
% % colormap(meshgrid(1:-1/63:0,1:3)'); 
% colormap(jet)
% 
% hold on;plot(Vphe_X(1:4:182),fe(1:4:182),'w.','Markersize',20);
% shading interp;
% title(['Radial rc=',num2str(last_fnu)])
% xlabel('Phase Velocity, v_{ph} (m/s)');ylabel('Frequency, f (Hz)');
% axis([150 350 10 100]);
% % title('Dispersion curve')
% set(gcf,'color','w');
% set(gca,'fontsize',18);

% %% Plot in fx domain
% 
% figure;
% subplot(121)
% colormap(jet)
% imagesc(x,f,abs(U_fZ/100),'Interpolation', 'bilinear')
% title('Vertical')
% xlabel('Offset (m)');ylabel('Frequency, f (Hz)');
% axis([10 50 5 80]);
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% colorbar
% subplot(122)
% colormap(jet)
% imagesc(x,f,abs(U_fX/100),'Interpolation', 'bilinear')
% title('Radial')
% xlabel('Offset (m)');ylabel('Frequency, f (Hz)');
% axis([10 50 5 80]);
% set(gca,'FontSize',18);
% colorbar
% set(gcf,'color','w');
% %% Plot in lx domain
% 
% figure,imagesc(x,le,abs(U))
% xlabel('Offset (m)');ylabel('Wavelength (m)');
% set(gca,'FontSize',30);
% set(gcf,'color','w');
 %% Perform f-l transform, f-l
 [f,l,A,U]=fl(yZ,dt,x,fmin,fmax,df,lmin,lmax,dl);
 [f,l,A2,U]=fl(yX,dt,x,fmin,fmax,df,lmin,lmax,dl); 

% %% Plot in fl domain
fe=f;
le=zeros(length(fe),1);
for i=1:length(fe);
   [Amax Ai]=max(A(i,:));
   le(i)=l(Ai);
    A(i,:) = A(i,:)/max(A(i,:));
    [A2max A2i]=max(A2(i,:));
   le(i)=l(Ai);
    A2(i,:) = A2(i,:)/max(A2(i,:));
end
% figure;
% imagesc(l,f,A);
% colormap(jet); 
% xlabel('Wavelength (m)');
% ylabel('Frequency, f (Hz)');
% axis([1 50 5 100]);
% hold on;
% plot(le,fe,'r.','Markersize',20);hold on;
% set(gca,'FontSize',30);
% set(gcf,'color','w');
% %% Plot in lv domain
% 
% figure; 
% plot(Vphe, le, 'b.','Markersize',20); hold on; 
% xlabel('Phase Velocity(m/s)'); ylabel('Wavelength (m)');  set(gca,'Ydir','reverse');
% axis([150 250 0 40]);
% legend('homo')
% set(gca,'FontSize',30);
% set(gcf,'color','w');



% %%  spectral amplitude normalization for each frequency
% [f,l,X_fl,R2_X]=fx_to_fl(U_fX,f,x',lmin,lmax,dl);
% [f,l,Z_fl,R2_Z]=fx_to_fl(U_fZ,f,x',lmin,lmax,dl);
% AX_fl=abs(X_fl); AZ_fl=abs(Z_fl);
% 
% for i=1:length(f);
%     [AXmax_fl AXi_fl]=max(AX_fl(i,:)); lXpeak(i)=l(AXi_fl);
%     [AZmax_fl AZi_fl]=max(AZ_fl(i,:)); lZpeak(i)=l(AZi_fl);
%     AXn_fl(i,:) = AX_fl(i,:)/AXmax_fl;
%     AZn_fl(i,:) = AZ_fl(i,:)/AZmax_fl; 

% end 
% %% Peform l-v transform, lv
%  [l,v,AAX]=lv(yX,dt,x,f,l,v);
%  [l,v,AAZ]=lv(yZ,dt,x,f,l,v);
% 
% 
% 
% figure;
% subplot(141);
% imagesc(l,f,AXn_fl);
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs; colormap(cmap); colorbar('southoutside');
% hold on;
% plot(lZpeak,f,'r.','Markersize',20);
% shading interp;
% title('Vertical','FontSize',30); xlabel('Wavelength(m)','FontSize',18); ylabel('frequency(Hz)','FontSize',30);set(gca,'FontSize',30);set(gcf,'color','w');
% % legend('b/h_{t}=1')
% legend('homo');
% subplot(142);
% imagesc(l,f,AZn_fl);
% colorbar('southoutside');
% hold on;
% plot(lXpeak,f,'r.','Markersize',20);
% title('Radial','FontSize',30); xlabel('Wavelength(m)','FontSize',18); ylabel('frequency(Hz)','FontSize',30);set(gca,'FontSize',30);set(gcf,'color','w');
% % legend('b/h_{t}=1')
% legend('homo');
% 
% figure;
% subplot(1,4,1)
% imagesc(v,l,AAZ);
% % cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs; colormap(cmap);
% % colorbar('southoutside');
% caxis([0 1]);
% hold on;
% plot(Vphe_Z,lZpeak,'w.','Markersize',20)
% shading interp;
% 
% title(['Vertical vc=',num2str(last_fnu)],'FontSize',24); xlabel('Phase Velocity, v_{ph} (m/s)','FontSize',18); ylabel('Wavelength(m)','FontSize',14);set(gca,'FontSize',18);set(gcf,'color','w');
% % legend('b/h_{t}=0.5');
% axis([100 350 0.5 10]);
% subplot(1,4,2)
% imagesc(v,l,AAX);
% % cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs; colormap(cmap);
% % colorbar('southoutside');
% caxis([0 1]);
% hold on;
% plot(Vphe_X,lXpeak,'w.','Markersize',20)
% title(['Radial rc=',num2str(last_fnu)],'FontSize',24); xlabel('Phase Velocity, v_{ph} (m/s)','FontSize',18); ylabel('Wavelength(m)','FontSize',14);set(gca,'FontSize',18);set(gcf,'color','w');
% % legend('b/h_{t}=0.5');
% axis([100 350 0.5 10]);
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
%% Plot ricker
% Create the wavelet and shift in time if needed
% figure;
% subplot(121)
% [rw,t] = ricker(30,N,0.00008);
% plot(t,rw)
% xlabel('Time (s)')
% ylabel('Amplitude')
% set(gca,'FontSize',24);
% set(gcf,'color','w');
% 
% subplot(122)
% RW=fft(rw);
% RW1=abs(RW)/max(abs(RW));
% RW2=RW1(1:182);
% plot(f,RW2)
% 
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% set(gca,'FontSize',24);
% set(gcf,'color','w');


% sRVSR, RVPD, RVSR

CWU=conj(U_fZ).*U_fX;
% RVSR = abs(U_fX)./abs(U_fZ); 
RVSR = abs(U_fZ)./abs(U_fX); 
RVPD=angle(CWU); 
% sRVSR=sign(imag(CWU)).*RVSR; 
sRVSR=sign(imag(CWU)).*RVSR; 
%% Subplot: sRVSR, RVPD, RVSR
figure;
subplot(133),pcolor(x,f,sRVSR); 
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',18,'Ydir','reverse');
set(gcf,'color','w');
% colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 
colormap(gca,jet);colorbar;caxis([-2 2]);
title('sRVSR'); xlabel('Offset (m)'); ylabel('Frequency, f (Hz)');

subplot(132),pcolor(x,f,RVPD);
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',18,'Ydir','reverse');
% colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 

colormap(gca,jet);colorbar;caxis([-pi pi]);
title('RVPD'); xlabel('Offset (m)'); ylabel('Frequency, f (Hz)');

subplot(131),pcolor(x,f,RVSR);
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',18,'Ydir','reverse','ColorScale','log');
% colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 

colormap(gca,gray);colorbar;caxis([-2 2]);
title('RVSR'); xlabel('Offset (m)'); ylabel('Frequency, f (Hz)');


% %% backscatter
% % -------------------------
% ybackscatterZ=extract_backscatter(yZ, dt); 
% figure; % Examine the original waveforms
% subplot(211); 
% wigb(yZ,1.,x,t);
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% ylabel('Time (s)'); xlabel('Offset(m)');
% title('Vertical Original Signal');
% subplot(212); 
% wigb(ybackscatterZ,1.,x,t);
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% ylabel('Time (s)'); xlabel('Offset (m)');
% title('Vertical extracted backscatter'); 
% 
% 
% ybackscatterX=extract_backscatter(yX, dt); 
% figure; % Examine the original waveforms
% subplot(211); 
% wigb(yX,1.,x,t);
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% ylabel('Time (s)'); xlabel('Offset (m)');
% title('Radial Original Signal');
% subplot(212); 
% wigb(ybackscatterX,1.,x,t);
% set(gca,'FontSize',18);
% set(gcf,'color','w');
% ylabel('Time (s)'); xlabel('Offset (m)');
% title('Radial extracted backscatter'); 
