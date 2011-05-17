export MYBASE=$HOME/pymol/movies
export LD_LIBRARY_PATH=${MYBASE}/bin:${LD_LIBRARY_PATH}
export PATH=${MYBASE}/bin:${PATH}
export PYM=pymol1.3
export FREEMOL=${MYBASE}/bin/freemol
export CLUSTER_USERNAME=abelew
export CLUSTER_HOSTNAME=login.deepthought.umd.edu
export CLUSTER_DIR=lustre
export CLUSTER_QSTAT=/usr/local/torque/bin/qstat
export MENCODER_OPTS="-ovc -lavc -lavcopts vcodec=mjpeg:vhq:psnr -noskip"