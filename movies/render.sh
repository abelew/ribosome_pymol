#!/bin/bash
cd ~/lustre
rsync -avz --progress --stats . abelew@login.deepthought.umd.edu:lustre/
ssh abelew@login.deepthought.umd.edu lustre/submit.sh

