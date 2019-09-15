#!/bin/bash

source /data/users/jabecass/dl_tools_centos/miniconda3_activate
source activate my_python36_centos

./src/get_feature_table.py /data/tmp/torque_logs/jabecass/ /data/tmp/torque6_logs/jabecass/