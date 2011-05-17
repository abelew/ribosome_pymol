#!/bin/bash
. render_options.bash
mencoder mf://*.png -o ../$DIRECTORY.avi $MENCODER_OPTS
