function [y2ifft]=extract_backscatter(y2,dt)
% Read seismic data
%[y, dt]=read_org('mc403_sx.org');
%load y.txt; dt=0.002; y=y(1:1000,:);

[N, Ch]=size(y2); % Note N and Ch must be even number here

Y2=fft(y2, N+1); Y2=Y2(1:N/2+1,:); 

% ------Apply spatial window---------
W=[]; 
w=kaiser(Ch,4)';
%w=ones(1,Ch);
for i=1:length(Y2)
   W=[W; w];
end
Y2=Y2.*W; 
% ----------------------------------

YY2=fft((Y2).',Ch+1).'; 

%% backward filter
%YY2(:,2:Ch/2+1)=zeros(N/2+1,Ch/2); % mute backward propagating waves
YY2(:,Ch/2+2:end)=zeros(N/2+1,Ch/2); % mute foreward propagating waves

% Return to time domain
Y2_filt=ifft(YY2.').'; Y2_filt=Y2_filt(:,1:Ch); 
Y2_filt=[Y2_filt; conj(flipud(Y2_filt(2:end,:)))]; 
y2ifft=ifft(Y2_filt);y2ifft=real(y2ifft(1:N,:));  
