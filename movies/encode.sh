#!/bin/bash
cd ~/lustre
rsync -avz --progress --stats abelew@login.deepthought.umd.edu:lustre/ .

echo "Please name the directory to encode."
read -e DIRECTORY
cd $DIRECTORY
mencoder mf://*.png -o ../$DIRECTORY.avi -ovc lavc -lavcopts vcodec=mjpeg:vhq:psnr -noskip

