
img_name='teddy_ear.bmp';

scribs_img_name='teddy_ear_m.bmp';

I = im2double(imread(img_name));
mI = im2double(imread(scribs_img_name));


levels_num=1;
active_levels_num=1;


sig=0.1^5;

alpha = runMattingfun(I,mI);