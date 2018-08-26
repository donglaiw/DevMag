% need libf2c and lapack blas
DIR_f2c = '../f2c/'; % folder containing f2c.h
DIR_clapack = '/opt/local/include/';
eval([' mex -g CFLAGS="\$CFLAGS -std=c99" -I' DIR_f2c ' -I' DIR_clapack ' -L' DIR_f2c ' elsd_mex_line.c valid_curve.c process_curve.c process_line.c write_svg.c -llapack -lblas  -lf2c -lm'])