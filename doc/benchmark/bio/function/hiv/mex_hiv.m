mex -output hiv_con_mex GCC='/home/maeda/local/bin/gcc4.7.4' COMPFLAGS='-Wall -O2 -I/home/maeda/local/include' source/Entrance.c source/SSR_hiv.c source/hiv.c source/cvSimulator.c /home/maeda/local/lib/libsundials_cvode.a /home/maeda/local/lib/libsundials_nvecserial.a