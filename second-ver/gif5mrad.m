clear all ; clc


cd D:\study\Hank\onelayercavitysize5mxyz20022525ht5mh2mmeshsize0.7m0.8sec0point5\
%// Image source: http:\\giantbomb.com
[A,map1] = rgb2ind(imread('radial lv 10~15m.jpg'),256);
[B,map2] = rgb2ind(imread('radial lv 12~17m.jpg'),256);
[C,map3] = rgb2ind(imread('radial lv 14~19m.jpg'),256);
[D,map4] = rgb2ind(imread('radial lv 16~21m.jpg'),256);
[E,map5] = rgb2ind(imread('radial lv 18~23m.jpg'),256);
[F,map6] = rgb2ind(imread('radial lv 20~25m.jpg'),256);
[G,map7] = rgb2ind(imread('radial lv 22~27m.jpg'),256);
[H,map8] = rgb2ind(imread('radial lv 24~29m.jpg'),256);
[I,map9] = rgb2ind(imread('radial lv 26~31m.jpg'),256);
[J,map10] = rgb2ind(imread('radial lv 28~33m.jpg'),256);
[K,map11] = rgb2ind(imread('radial lv 30~35m.jpg'),256);
[L,map12] = rgb2ind(imread('radial lv 32~37m.jpg'),256);
[M,map13] = rgb2ind(imread('radial lv 34~39m.jpg'),256);
[N,map14] = rgb2ind(imread('radial lv 36~41m.jpg'),256);
[O,map15] = rgb2ind(imread('radial lv 38~43m.jpg'),256);
[P,map16] = rgb2ind(imread('radial lv 40~45m.jpg'),256);
[Q,map17] = rgb2ind(imread('radial lv 42~47m.jpg'),256);
[R,map18] = rgb2ind(imread('radial lv 44~49m.jpg'),256);
[S,map19] = rgb2ind(imread('radial lv 45~50m.jpg'),256);
ImageCell = {A;B;C;D;E;F;G;H;I;J;K;L;M;N;O;P;Q;R;S};
%// Just to show what the images look like (I removed spots to make sure there was an animation created):
%// Create file name.
FileName = 'radial lv.gif';
for k = 1:numel(ImageCell)
    
    if k ==1
        
        %// For 1st image, start the 'LoopCount'.
        imwrite(ImageCell{k},map1,FileName,'gif','LoopCount',Inf,'DelayTime',0.5);
    elseif k ==2
        imwrite(ImageCell{k},map2,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==3
        imwrite(ImageCell{k},map3,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==4
        imwrite(ImageCell{k},map4,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==5
        imwrite(ImageCell{k},map5,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==6
        imwrite(ImageCell{k},map6,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==7
        imwrite(ImageCell{k},map7,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==8
        imwrite(ImageCell{k},map8,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==9
        imwrite(ImageCell{k},map9,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==10
        imwrite(ImageCell{k},map10,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==11
        imwrite(ImageCell{k},map11,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==12
        imwrite(ImageCell{k},map12,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==13
        imwrite(ImageCell{k},map13,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==14
        imwrite(ImageCell{k},map14,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==15
        imwrite(ImageCell{k},map15,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==16
        imwrite(ImageCell{k},map16,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==17
        imwrite(ImageCell{k},map17,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==18
        imwrite(ImageCell{k},map18,FileName,'gif','WriteMode','append','DelayTime',0.5);
    elseif k ==19
        imwrite(ImageCell{k},map19,FileName,'gif','WriteMode','append','DelayTime',0.5);
    end
    
end