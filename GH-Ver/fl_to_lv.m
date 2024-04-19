function [l,v, AA]=fl_to_lv(A, f, l, v) 
% This function transform fl to lv domain
% 
%Latest updated by CP Lin 2021/12/30

% Spectral amplitdue normalization for each frequency
for i=1:length(f);
    [Amax Ai]=max(A(i,:)); 
    A(i,:) = A(i,:)/Amax;
end

%% Genetrating lamda_velocity spectrum from fl spectrum by interpolation
[ll ff]=meshgrid(l,f); 
ln = reshape(ll,[],1); 
fn = reshape(ff, [], 1); 
An=reshape(A, [],1); 
vn=ln.*fn; 
[vv,llz]=meshgrid(v,l); 
AA=griddata(vn, ln, An, vv, llz, 'cubic'); % among 'linear', 'nearest', 'natural', 'cubic', or 'v4'.