clear all,close all,clc
cd('D:\模擬\學姊\');
% In the practicing now:
% Z: vertical direction, 
% X: radial direction, 
% Y: transversal direction

%% Simulated files combine (資料讀入)
first_fnu = 1; last_fnu = 20; % The start and end station for analysis
NF=last_fnu-first_fnu+1; % Number of files


fn1='DB.X';
Direction_X='X';fn3_X=['.FX',Direction_X,'.semd'];
Direction_Y='Y';fn3_Y=['.FX',Direction_Y,'.semd'];
Direction_Z='Z';fn3_Z=['.FX',Direction_Z,'.semd'];

for i=1:1:NF % 
    fn2=num2str(first_fnu+i-1); % The following lines putting the pieces of file names together
    FileName_X=[fn1 fn2 fn3_X]; % 彙整檔案名 (X)
    FileName_Y=[fn1 fn2 fn3_Y]; % (Y)
    FileName_Z=[fn1 fn2 fn3_Z]; % (Z)

    RawData_X=load(FileName_X); % 讀入X_Y_Z檔案
    RawData_Y=load(FileName_Y);
    RawData_Z=load(FileName_Z);
    
    % 不太懂 是為了讓資料不要是負的嗎? 為何需要這步
    y_X(:,i)=-RawData_X(:,2); % Add negitive sign for direction correction
    y_X_norm(:,i) = y_X(:,i)/max(abs(y_X(:,i))); % 正規化
    y_Y(:,i)=-RawData_Y(:,2); % Add negitive sign for direction correction
    y_Y_norm(:,i) = y_Y(:,i)/max(abs(y_Y(:,i))); % 正規化
    y_Z(:,i)=RawData_Z(:,2);
    y_Z_norm(:,i) = y_Z(:,i)/max(abs(y_Z(:,i)));% 正規化

end

[N, Ch]=size(y_Z); % 獲得X、Y軸數量 X:N、Y:Ch
% The offset for analysis (Carefull)
dx=0.1; % Geophone spacing
x0=0.5; % distance between source & 1st receiver
x=[x0:dx:x0+(Ch-1)*dx]'; % 標出每個測站位置

t = RawData_X(:,1); dt=t(2)-t(1); % time step t:每筆資料的時間 dt:時間間隔 
% [m,t0_ind] = min(abs(t)); % Time correction, the simulation time start from negitive time 
% t = t(t0_ind:end); N = length(t);
% y_X = y_X(t0_ind:end,:); y_Y = y_Y(t0_ind:end,:); y_Z = y_Z(t0_ind:end,:);
N = length(t); t=[0:dt:(N-1)*dt]; % [起始:跌代值:終值]

%% Figure plot in time domain
figure,
subplot(131),wigb(y_Z_norm,1,x,t); % plot seismic data in wiggle format
set(gca,'FontSize',10,'Ydir','reverse');ylabel('Time (s)'); xlabel('Offset (m)'); title('Vertical component')
subplot(132),wigb(y_X_norm,1,x,t);
set(gca,'FontSize',10,'Ydir','reverse');ylabel('Time (s)'); xlabel('Offset (m)'); title('Radial component')
subplot(133),wigb(y_Y_norm,1,x,t);
set(gca,'FontSize',10,'Ydir','reverse');ylabel('Time (s)'); xlabel('Offset (m)'); title('Transversal component')

%% Zero padding to 1 sec
% 目前這邊的理解是 將資料拓展成以頻率為基準的資料集
N = round(1/dt);t = [0:dt:(N-1)*dt]; % 將頻率算出 N:f
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
fmin=0; fmax=14000; df=10;
% Analysis velocity range
vmin = 0; vmax=2000; dv=10;
% Analysis wave length range
lmin=0; lmax=100; dl=0.5;
% Analysis wave number range
kmin=-80; kmax=80; dk=0.05;

% Peform f-v transform, f-v

[f,v,A_Z,Z_fx]=fv(y_Z,dt,x,fmin,fmax,df,vmin,vmax,dv);
[f,v,A_X,X_fx]=fv(y_X,dt,x,fmin,fmax,df,vmin,vmax,dv);
[f,v,A_Y,Y_fx]=fv(y_Y,dt,x,fmin,fmax,df,vmin,vmax,dv);

% Peform f-k transform, f-k
% [f,v,A_Z,Z_fx]=fk(y_Z,dt,x,fmin,fmax,df,kmin,kmax,dk);
% [f,v,A_X,X_fx]=fk(y_X,dt,x,fmin,fmax,df,kmin,kmax,dk);
% [f,v,A_Y,Y_fx]=fk(y_Y,dt,x,fmin,fmax,df,kmin,kmax,dk);

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

% for i=1:length(fe);
%    [Amax Ai]=max(Ak_Z(i,:));
%    Vphe_Z(i)=v(Ai);
%    An_Z_k(i,:) = Ak_Z(i,:)/Amax;
% 
%    [Amax Ai]=max(Ak_X(i,:));
%    Vphe_X(i)=v(Ai);
%    An_X_k(i,:) = Ak_X(i,:)/Amax;
% 
%    [Amax Ai]=max(Ak_Y(i,:));
%    Vphe_Y(i)=v(Ai);
%    An_Y_k(i,:) = Ak_Y(i,:)/Amax;
% end


%% Radial to vertical relation (sRVSR, RVPD and RVSR)
CWU=conj(Z_fx).*X_fx;
RVSR = abs(X_fx)./abs(Z_fx); 
RVPD=angle(CWU); 
sRVSR=sign(imag(CWU)).*RVSR;  

% Subplot: sRVSR, RVPD, RVSR
figure;
% subplot(131),
pcolor(x,f,sRVSR); 
shading interp;
set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse');
colormap(gca,jet);colorbar;caxis([-2 2]);
title('sRVSR'); xlabel('x (m)'); ylabel('f (Hz)');

% subplot(132),pcolor(x,f,RVPD);
% shading interp;
% set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse');
% colormap(gca,jet);colorbar;caxis([-pi pi]);
% title('RVPD'); xlabel('x (m)'); ylabel('f (Hz)');
% 
% subplot(133),pcolor(x,f,RVSR);
% shading interp;
% set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
% colormap(gca,gray);colorbar;caxis([-2 2]);
% title('RVSR'); xlabel('x (m)'); ylabel('f (Hz)');

%% Plot in f-x domain
% figure,
% subplot(131),pcolor(x,f,abs(Z_fx));
% shading interp;
% set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
% % colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 
% colormap(gca,jet);colorbar;%caxis([-2 2]);
% title('Vertical'); xlabel('x (m)'); ylabel('f (Hz)');
% 
% subplot(132),pcolor(x,f,abs(X_fx));
% shading interp;
% set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
% % colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 
% colormap(gca,jet);colorbar;%caxis([-2 2]);
% title('Radial'); xlabel('x (m)'); ylabel('f (Hz)');
% 
% subplot(133),pcolor(x,f,abs(Y_fx));
% shading interp;
% set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse','ColorScale','log');
% % colormap(gca, cbrewer2('RdGy'));colorbar; caxis([-2 2]); 
% colormap(gca,jet);colorbar;%caxis([-2 2]);
% title('Transversal'); xlabel('x (m)'); ylabel('f (Hz)');


%% Plot in fv domain
figure;
% subplot(131),imagesc(f,v,An_Z);colormap(jet);
imagesc(f, v, An_Z);colormap(jet); % 先別加 subplot，圖會變小張 4/13
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);
hold on;plot(Vphe_Z,fe,'w.','Markersize',10);
ylabel('Phase Velocity, V_{ph} (m/s)');xlabel('Frequency, f (Hz)'); title('Vertical')
set(gca,'fontsize',14,'YDir', 'normal');

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
% figure;
% subplot(131),imagesc(v,l,A_lv_Z);colormap(jet);
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);caxis([0 1]);
% hold on;plot(Vphe_Z,Vphe_Z./fe,'w.','Markersize',10);
% xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Wavelength, \lambda (m)'); title('Vertical')
% set(gca,'fontsize',14);
% 
% subplot(132),imagesc(v,l,A_lv_X);colormap(jet);
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);caxis([0 1]);
% hold on;plot(Vphe_X,Vphe_X./fe,'w.','Markersize',10);
% xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Wavelength, \lambda (m)'); title('Radial')
% set(gca,'fontsize',14);
% 
% subplot(133),imagesc(v,l,A_lv_Y);colormap(jet);
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);caxis([0 1]);
% hold on;plot(Vphe_Y,Vphe_Y./fe,'w.','Markersize',10);
% xlabel('Phase Velocity, V_{ph} (m/s)');ylabel('Wavelength, \lambda (m)'); title('Transversal')
% set(gca,'fontsize',14);


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
