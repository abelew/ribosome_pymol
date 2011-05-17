#!/bin/bash
. render_options.bash

ssh ${CLUSTER_USERNAME}@${CLUSTER_HOSTNAME} $CLUSTER_QSTAT
