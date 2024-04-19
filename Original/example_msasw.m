% This is an example of MSASW analysis

% Read seismic data
[Data, dt]=read_org('Mc403_sx.org');

% Perform mutichanel SASW analysis
df=1.0;fmin=5; fmax=80; % Frequency Range
dx=1; % Geophone spacing
[N, Ch]=size(Data);x1=10; x=[x1:dx:x1+(Ch-1)*dx]';

[fe, Vphe, R2, UangleU, Uangle]=MSASW(Data,dt,x,fmin,fmax,df);


% Plotting the f-x phase image and dispersion curve
figure;
imagesc(x,fe,UangleU);
%colormap(meshgrid(1:-1/63:0,1:3)');
c1=0; c2=0; c3=0; cs=1;
cmap=[1:-(1-c1)/63:c1; 1:-(1-c2)/63:c2; 1:-(1-c3)/63:c3]';
cmap=cmap.^cs;
colormap(cmap);
colorbar;
xlabel('x, (m)');
ylabel('Frequency, f (Hz)');
hold on;

% Plots
figure
plot(fe,Vphe,'b*');hold on; 
set(gca,'FontSize',16,'LineWidth',2);
ylabel('V_p_h (m/s)');
xlabel('f (Hz)');
axis([0 80 100 400]);
