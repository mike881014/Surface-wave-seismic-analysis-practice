% This is an example of SASW analysis

% Read seismic data
[Data, dt]=read_org('Mc403_sx.org');

% Perform SASW analysis
%[f, Vph]=sasw(Data(:,5), Data(:,13), dt, 8,1/3,2);
%[f1, Vph1]=sasw(Data(:,13), Data(:,21), dt, 8,1/3,2);
%[f2, Vph2]=sasw(Data(:,15), Data(:,23), dt, 8,1/3,2);

[f, Vph]=sasw(Data(:,13), Data(:,21), dt, 8,1/3,1.5);
[f1, Vph1]=sasw(Data(:,8), Data(:,24), dt, 16,1/3,1.5);
[f2, Vph2]=sasw(Data(:,2), Data(:,22), dt, 20,1/3,1.5);

% Plots
figure
plot(f,Vph,'k*');hold on; 
plot(f1,Vph1,'b.'); 
plot(f2,Vph2, 'ro');
set(gca,'FontSize',16,'LineWidth',2);
xlabel('f (Hz)', 'FontSize', 16);
ylabel('Vs (m/s)', 'FontSize', 16); 
legend('\DeltaX=8', '\DeltaX=16', '\DeltaX=24');
axis([15 80 100 400]);