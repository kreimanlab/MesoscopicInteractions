#!/bin/bash

# Figure W2-*A
matlab -nodesktop -nosplash -nosoftwareopengl -r "figure_t14d1_raw; exit();"

# Figure W2-*B
matlab -nodesktop -nosplash -nosoftwareopengl -r "figure_t14d1; exit();"

# Figure W2-*C
matlab -nodesktop -nosplash -nosoftwareopengl -r "figure_t14d1_150; exit();"

# Figure S?
#matlab -nodesktop -nosplash -nosoftwareopengl -r "figure_t14_150; exit();"

# Figure W3
matlab -nodesktop -nosplash -nosoftwareopengl -r "figure_t14_allatl; exit();"
