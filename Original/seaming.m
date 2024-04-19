function [y,dt]=seaming(first_fnu,last_fnu)

% [y, dt]=seaming(first_fnu,last_fnu)
% This function seaming the phase difference at overlapped offsets
% 		y = corrected data y(t,x)
% 		first_fnu= first file number of the common receiver data
%		last_fnu= last file number of the common receiver data
% Written by C-P Lin, Last modified 11/14/2003.

% Read OYO data mc???_sx.org
NF=last_fnu-first_fnu+1;
   fn1='mc';
	fn3='_sx.org';
for i=1:NF
   fn2=num2str(first_fnu+i-1);
   fn2_length=length(fn2);
   if fn2_length<3
      fn2=[num2str(zeros(1,3-fn2_length),'%1d') fn2];
   end
   FileName=[fn1 fn2 fn3];
   [y(:,:,i), dt]=read_org(FileName);
end
t=[0:dt:(length(y(:,1,1))-1)*dt]';


% Seaming in frequency domain
[M N NF]=size(y);
Y=fft(y);
Yangle=angle(Y);

% correct Yangle, Y, and y
j=sqrt(-1);
for i=1:NF-1
   dangle=Yangle(:,N,i)-Yangle(:,1,i+1);
   dangle=repmat(dangle,1,N);
   shift=exp(j*dangle);
   Y(:,:,i+1)=Y(:,:,i+1).*shift;
   Yangle(:,:,i+1)=angle(Y(:,:,i+1));
end
Y=[reshape(Y(:,1:N-1,:),M,(N-1)*NF) Y(:,N,NF)];
y=real((ifft(Y)));
[M Ch]=size(Y);



