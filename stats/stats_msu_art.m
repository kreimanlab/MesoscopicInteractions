close all;
clear;

a = '/media/jerry/internal/data/h5_notch20/art_nosz/mSu_art.h5';
art_idx = h5read(a,'/art_idx');



imagesc(art_idx);




%After LOC was established, neural activity was recorded for ~25 min for the ketamine and medetomidine-induced anesthesia experiment and ~10 min for the propofol-induced anesthesia experiment. Heart rate and breathing were monitored carefully throughout the length of the experiment. 