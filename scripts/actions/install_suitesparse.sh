wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v5.9.0.tar.gz 
tar -zxvf v5.9.0.tar.gz && cd SuiteSparse-5.9.0
make library -j4
sudo make install INSTALL=/usr/local/suitesparse

