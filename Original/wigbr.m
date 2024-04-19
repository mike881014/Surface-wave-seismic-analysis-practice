function wig (Data,scal,x,t,amx)
%WIGB    WIGB(Data,scal,x,t,amx) to plot seismic data in wiggle format.
%  'Data' is input matrix data(t,x)
%	 x and t are vectors with offset and time.
%	 If only 'Data' is enter, 'scal,x,t,amn,amx' are decided automatically; 
%	 otherwise, 'scal' is a scalar; 'x, t' are vectors for annotation in 
%	offset and time, amx are the amplitude range.


if nargin == 0, nx=10;nz=10; Data = rand(nz,nx)-0.5; end;

[nz,nx]=size(Data);

trmx= max(abs(Data));
if (nargin <= 4); amx=mean(trmx);  end;
if (nargin <= 2); x=[1:nx]; t=[1:nz]; end;
if (nargin <= 1); scal =1; end;

if nx <= 1; disp(' ERR:PlotWig: nx has to be more than 1');return;end;

 % take the average as dx
	dx1 = abs(x(2:nx)-x(1:nx-1));
 	dx = median(dx1);

 dz=t(2)-t(1);
 xmx=max(max(Data)); xmn=min(min(Data)); 

 if scal == 0; scal=1; end;
 Data = Data * dx /amx; 
 Data = Data * scal;

 %fprintf(' PlotWig: data range [%f, %f], plotted max %f \n',xmn,xmx,amx);
 
% set display range 
x1=min(x)-2.0*dx; x2=max(x)+2.0*dx;
z1=min(t)-dz; z2=max(t)+dz;
 
set(gca,'NextPlot','add','Box','on', ...
  'XLim', [z1 z2], ...
  'YDir','reverse', ...
  'YLim',[x1 x2]);
 

	fillcolor = [0 0 0];
	linecolor = [0 0 0];
	linewidth = 0.1;

	t=t'; 	% input as row vector
	zstart=t(1);
   zend  =t(nz);
   
for i=1:nx,
   
  if trmx(i) ~= 0;    % skip the zero traces
	tr=Data(:,i); 	% --- one scale for all section
  	s = sign(tr) ;
  	i1= find( s(1:nz-1) ~= s(2:nz) );	% zero crossing points
	npos = length(i1);


	%12/7/97 
	zadd = i1 + tr(i1) ./ (tr(i1) - tr(i1+1)); %locations with 0 amplitudes
	aadd = zeros(size(zadd));

	[zpos,vpos] = find(tr >0);
	[zz,iz] = sort([zpos; zadd]); 	% indices of zero point plus positives
	aa = [tr(zpos); aadd];
	aa = aa(iz);

	% be careful at the ends
		if tr(1)>0, 	a0=0; z0=1.00;
		else, 		a0=0; z0=zadd(1);
		end;
		if tr(nz)>0, 	a1=0; z1=nz; 
		else, 		a1=0; z1=max(zadd);
		end;
			
	zz = [z0; zz; z1; z0];
 	aa = [a0; aa; a1; a0];
		

	zzz = zstart + zz*dz -dz;
   
   patch(zzz, aa+x(i), fillcolor);

	line( 'Color',[1 1 1],'EraseMode','background',  ...
         'Xdata', [zstart zend], 'Ydata',x(i)+[0 0]); % remove zero line

      
%'LineWidth',linewidth, ...
%12/7/97 'Xdata', x(i)+[0 0], 'Ydata',[z0 z1]*dz);	% remove zero line

	line( 'Color',linecolor,'EraseMode','background',  ...
	 'LineWidth',linewidth, ...
	 'Xdata',t , 'Ydata',tr+x(i));	% negatives line

   else % zeros trace
	line( 'Color',linecolor,'EraseMode','background',  ...
	 'LineWidth',linewidth, ...
         'Xdata',[zstart zend] , 'Ydata',[x(i) x(i)]);
   end;
end;





