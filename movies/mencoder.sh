#!/bin/sh
mencoder mf://*.png -ovc lavc -lavcopts vcodec=mjpeg:vhq:psnr -noskip -o mjpeg.avi
