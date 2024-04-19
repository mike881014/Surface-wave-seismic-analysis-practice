% This is a program for on-site quality control (Checking data)

% Read seismic data
[y, dt]=read_org('Mc403_SX.org');
dx=1; % Geophone spacing
x1=15; 

%y=fliplr(y);
%y=y(:,4:24);

[N, Ch]=size(y);
x=[x1:dx:x1+(Ch-1)*dx]';
t=0:dt:(N-1)*dt;

%---------------------------------------------------------------------------------------
% Peform f-p transform
fmin=1; fmax=100; df=1; % Frequency Range
pmin=1/500; pmax=1/70; dp=(pmax-pmin)/200; % Wavenumber Range
[f,p,A,U]=fp(y,dt,x,fmin,fmax,df,pmin,pmax,dp);

% Find the dispersion curve from the overone image
pe=zeros(length(f),1);

for i=1:length(f);
   [Amax Ai]=max(A(i,:));
   pe(i)=p(Ai);
end

%---------------------------------------------------------------------------------------
% Perform MSASW regression of phi-x
p_unwrap = 1.0;
UangleU=unwrap(angle(U),[p_unwrap*pi],2); % Unwrap phase angle

pe2=zeros(length(f),1);
R2=zeros(length(f),1);
for i=1:length(f)
   phi=UangleU(i,:);
   a=1; b=Ch; % Effective channel from a to b
   %b=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis
   [b bint r rint stats]=regress(x(a:b),[phi(a:b)' ones(b-a+1,1)]); % Regression analysis with R2 statistic
   %Vphe2(i)=abs(b(1))*2*pi*f(i);
   pe2(i)=1./(2*pi*f(i)*abs(b(1)));
   R2(i)=stats(1); % R2 statistic
end

%---------------------------------------------------------------------------------------
% Plot seismic data in t-x domain
figure
subplot(141)
wigb(y,1.,x,t);
set(gca,'FontSize',10);
ylabel('Time (s)')
xlabel('Offset (m)')
set(gca,'Ydir','reverse');

%plot in f-x domain amaplitude
subplot(142);
wigb(abs(U),1.,x,f);
set(gca,'FontSize',10);
ylabel('frequency (Hz)')
xlabel('Offset (m)')
set(gca,'Ydir','reverse');

% plot in f-x domain R2 of phi-x regression
subplot(143);
plot(R2,f,'b'); hold on;
axis([0 1.0 0 100]);
set(gca,'Ydir','reverse','FontSize',10);
ylabel('Frequency (Hz)','FontSize',10);
xlabel('R^2','FontSize',10);

% Plotting the f-p image and dispersion curve
subplot(144)
imagesc(p,f,A);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=1;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
xlabel('Slowness, p (s/m)');
ylabel('Frequency, f (Hz)');
hold on;
plot(pe,f,'r.');hold on;
plot(pe2,f,'bo');
legend('MWTSW','MSASW');
