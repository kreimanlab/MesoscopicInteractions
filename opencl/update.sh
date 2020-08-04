#!/bin/bash

rsync --exclude h5 --exclude remote --exclude xsub -avP $OPENCL_DIR/ .
