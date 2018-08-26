img_name='KNhair.bmp';
scribs_img_name='KNhair_m.bmp';

I = im2double(imread(img_name));
mI = im2double(imread(scribs_img_name));

alpha = runMattingfun(I,mI);

%runMatting