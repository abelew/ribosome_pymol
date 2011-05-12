#!/bin/bash
echo "This script will submit a ribosome movie session for encoding"
echo "It requires a couple environment variables to be set."
echo "1.  SESSIONDIR  :  a directory in which the pymol session file should live."
echo "2.  MOVIE_SCRIPT : A python script used to direct the movie."
#export MYBASE=/a/f20-fs1/data/dt-vol6/abelew
export MYBASE=~/lustre
export LD_LIBRARY_PATH=${MYBASE}/bin:${LD_LIBRARY_PATH}
export PATH=${MYBASE}/bin:$PATH
export PYM=${MYBASE}/bin/pymol1.4
echo "Type the name of the directory with your pymol session here."
echo "It should contain a single file named 'session.pse' inside it."
read -e NAME
export SESSIONNAME=$NAME
export SESSIONDIR=$MYBASE/$SESSIONNAME
export SESSION=${SESSIONDIR}/session.pse

if [ ! -f "$SESSION" ]; then
  echo "The file $SESSION does not exist."
  exit 1
fi

if [ "$FRAMES" = "" ]; then
  echo "The default number of frames in the movie is 480."
  echo "Change this by export FRAMES=###"
  export FRAMES=480
fi

if [ "$RAY" = "" ]; then
  export RAY=0
else
  export RAY=1
fi

if [ "$MOVIE_SCRIPT" = "" ]; then
  echo "The default movie script is render.py"
  echo "Change this by export MOVIE_SCRIPT=\"newscript.py\""
  export MOVIE_SCRIPT=$MYBASE/bin/local_render.py
else
  export MOVIE_SCRIPT=$MYBASE/bin/$MOVIE_SCRIPT
fi

cd $SESSIONDIR && $PYM -e -r $MOVIE_SCRIPT
