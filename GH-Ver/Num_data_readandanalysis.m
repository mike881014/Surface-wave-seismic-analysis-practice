clear all,close all,clc
cd('D:\模擬\Surface-wave-seismic-analysis-practice\GH-Ver');
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
    FileName_Z=[fn1 fn2 fn3_Z]; % (Z)

    RawData_Z=load(FileName_Z);
    
    % 不太懂 是為了讓資料不要是負的嗎? 為何需要這步
    y_Z(:,i)=RawData_Z(:,2);
    y_Z_norm(:,i) = y_Z(:,i)/max(abs(y_Z(:,i)));% 正規化

end

[N, Ch]=size(y_Z); % 獲得X、Y軸數量 X:N、Y:Ch
% The offset for analysis (Carefull)
dx=0.1; % Geophone spacing
x0=0.5; % distance between source & 1st receiver
x=[x0:dx:x0+(Ch-1)*dx]'; % 標出每個測站位置

t = RawData_Z(:,1); dt=t(2)-t(1); % time step t:每筆資料的時間 dt:時間間隔 
% [m,t0_ind] = min(abs(t)); % Time correction, the simulation time start from negitive time 
% t = t(t0_ind:end); N = length(t);
% y_X = y_X(t0_ind:end,:); y_Y = y_Y(t0_ind:end,:); y_Z = y_Z(t0_ind:end,:);
N = length(t); t=[0:dt:(N-1)*dt]; % [起始:跌代值:終值]

%% Figure plot in time domain
figure,
wigb(y_Z_norm,1,x,t); % plot seismic data in wiggle format
set(gca,'FontSize',10,'Ydir','reverse');ylabel('Time (s)'); xlabel('Offset (m)'); title('Vertical component')

%% Zero padding to 1 sec
% 目前這邊的理解是 將資料拓展成以頻率為基準的資料集
N = round(1/dt);t = [0:dt:(N-1)*dt]; % 將頻率算出 N:f
data_length = length(y_Z);
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
fmin=0; fmax=14000; df=1;
% Analysis velocity range
vmin = 0; vmax=2000; dv=1;
% Analysis wave length range
lmin=0; lmax=100; dl=1;
% Analysis wave number range
kmin=-80; kmax=80; dk=1;

% Peform f-v transform, f-v
[f,v,A_Z_v,Z_fx]=fv(y_Z,dt,x,fmin,fmax,df,vmin,vmax,dv);

% Peform f-k transform, f-k
[f,k,A_Z_k,Z_fk]=fk(y_Z,dt,x,fmin,fmax,df,kmin,kmax,dk);

% Perform f-l transform and then transform to l-v
[f,l,Al_Z]=fl(y_Z,dt,x,fmin,fmax,df,lmin,lmax,dl);[l,v, A_lv_Z]=fl_to_lv(Al_Z, f, l, v);


fe=f;
Vphe_Z_v=zeros(length(fe),1);
for i=1:length(fe);
   [Amax_v Ai_v]=max(A_Z_v(i,:));
   Vphe_Z_v(i)=v(Ai_v);
   An_Z_v(i,:) = A_Z_v(i,:)/Amax_v;
end

Vphe_Z_k=zeros(length(fe),1);
for i=1:length(fe);
   [Amax_k Ai_k]=max(A_Z_k(i,:));
   Vphe_Z_k(i)=k(Ai_k);
   An_Z_k(i,:) = A_Z_k(i,:)/Amax_k;
end


%% Radial to vertical relation (sRVSR, RVPD and RVSR) f-v version.
% CWU=conj(Z_fx).*X_fx;
% RVSR = abs(X_fx)./abs(Z_fx); 
% RVPD=angle(CWU); 
% sRVSR=sign(imag(CWU)).*RVSR;  
% 
% % Subplot: sRVSR, RVPD, RVSR
% figure;
% % subplot(131),
% pcolor(x,f,sRVSR); 
% shading interp;
% set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse');
% colormap(gca,jet);colorbar;caxis([-2 2]);
% title('sRVSR'); xlabel('x (m)'); ylabel('f (Hz)');

%% Radial to vertical relation (sRVSR, RVPD and RVSR) f-k version.
% CWU=conj(Z_fk).*X_fx;
% RVSR = abs(X_fk)./abs(Z_fk); 
% RVPD=angle(CWU); 
% sRVSR=sign(imag(CWU)).*RVSR;  
% 
% % Subplot: sRVSR, RVPD, RVSR
% figure;
% % subplot(131),
% pcolor(x,f,sRVSR); 
% shading interp;
% set(gca,'layer','top','Box','on', 'FontSize',12,'Ydir','reverse');
% colormap(gca,jet);colorbar;caxis([-2 2]);
% title('sRVSR'); xlabel('x (m)'); ylabel('f (Hz)');

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


%% Plot in f-v domain
figure;
% subplot(131),imagesc(f,v,An_Z);colormap(jet);
imagesc(f, v, An_Z_v);colormap(jet); % 先別加 subplot，圖會變小張 4/13
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);
hold on;plot(Vphe_Z_v,fe,'w.','Markersize',10);
ylabel('Phase Velocity, V_{ph} (m/s)');xlabel('Frequency, f (Hz)'); title('F-V Vertical')
set(gca,'fontsize',14,'YDir', 'normal');

%% Plot in f-k domain
figure;
% subplot(131),imagesc(f,v,An_Z);colormap(jet);
imagesc(f, k, An_Z_k);colormap(jet); % 先別加 subplot，圖會變小張 4/13
% cs=0.5; cmap=(meshgrid(1:-1/63:0,1:3)').^cs;colormap(cmap);
hold on;plot(Vphe_Z_k,fe,'w.','Markersize',10);
ylabel('Phase Velocity, V_{ph} (m/s)');xlabel('Frequency, f (Hz)'); title('F-K Vertical')
set(gca,'fontsize',14,'YDir', 'normal');

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
