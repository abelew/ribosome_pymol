#!/bin/bash
. render_options.bash
cd ${MYBASE}
rsync -avz --progress --stats ${CLUSTER_USERNAME}@${CLUSTER_HOSTNAME}:${CLUSTER_DIR}/ .

echo "Please name the directory to encode."
read -e DIR
export DIRECTORY=$DIR
cd $DIRECTORY
../bin/encode.sh