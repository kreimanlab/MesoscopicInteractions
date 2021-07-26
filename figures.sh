#!/bin/bash
mkdir stats/figures
echo source downloads/pyenv/bin/activate
source downloads/pyenv/bin/activate

export SUBJECTS_DIR='../data/coregistration'

currDir=`pwd`
cd stats

echo ""
echo "Figures"
echo ""

PS3=": "

select opt in "fig 1" "fig 2" "fig 3" "fig 4" "fig 5" "fig 6" "fig s1" "fig s2" "fig s3" "fig s4" "fig s5" "fig s6" "fig s7" "fig s8" "fig s9" "fig s10" "fig s11" "fig s12" "fig s13" "fig s14" quit; do

  case $opt in
    "fig 1")
      matlab -nodesktop -nosplash -r "figure_48surf;exit()"
      echo "Output: stats/figures/48surf"
      break
      ;;
    "fig 2")
      matlab -nodesktop -nosplash -r "figure_t1_2_3_4_6_9;exit()"
      echo "Output (A): stats/figures/T1"
      echo "Output (B): stats/figures/T9d1"
      echo "Output (C): stats/figures/T9"
      break
      ;;
    "fig 3")
      matlab -nodesktop -nosplash -r "figure_t8d1;figure_t7d1;exit()"
      echo "Output (A-B): stats/figures/T8d1"
      echo "Output (C-E): stats/figures/caret"
      echo "Output (F-I): stats/figures/T7d1"
      break
      ;;
    "fig 4")
      matlab -nodesktop -nosplash -r "figure_t19;figure_t16_dk;exit()"
      echo "Output: stats/figures/T16"
      break
      ;;
    "fig 5")
      matlab -nodesktop -nosplash -r "figure_t14_allatl;exit()"
      echo "Output: stats/figures/T14"
      break
      ;;
    "fig 6")
      matlab -nodesktop -nosplash -r "fig_cluster2_sworld_circle;figure_t22;exit()"
      echo "Output: stats/figures/fig_cluster2_sworld_circle_metric-1"
      echo "Output: stats/figures/T22"
      break
      ;;
    "fig 7")
      matlab -nodesktop -nosplash -r "figure_T21;exit()"
      break
      ;;
    "fig s1")
      mkdir figures/T4
      matlab -nodesktop -nosplash -r "figure_t4;exit()"
      echo "Output: stats/figures/T4"
      matlab -nodesktop -nosplash -r "dist_thresh;exit()"
      echo "Output: stats/figures/dist_thresh"
      break
      ;;
    "fig s2") #7
      break
      ;;
    "fig s3")
      matlab -nodesktop -nosplash -r "figure_t1_2_3_4_6_9;exit()"
      echo "Output (A): stats/figures/T9d1"
      echo "Output (B): stats/figures/T9"
      break
      ;;
    "fig s4")
      matlab -nodesktop -nosplash -r "figure_t1_2_3_4_6_9;exit()"
      echo "Output (A-D): stats/figures/T1d1"
      echo "Output (E): stats/figures/T9d1"
      echo "Output (F): stats/figures/T9"
      break
      ;;
    "fig s5")
      matlab -nodesktop -nosplash -r "figure_t1_2_3_4_6_9;exit()"
      echo "Output (A): stats/figures/T1d3"
      break
      ;;
    "fig s7")
      matlab -nodesktop -nosplash -r "figure_t1_2_3_4_6_9;exit()"
      echo "Output (A-D): stats/figures/T1d2/sub3_35_84_t1d2_one"
      break
      ;;
    "fig s8") #13
      matlab -nodesktop -nosplash -r "figure_t11_new;exit()"
      echo "Output (A-C): stats/figures/T11_new"
      #matlab -nodesktop -nosplash -r "s2fsaverage;exit()" # ielvis
      #requires MRI
      #echo "Output (D): data/coregistration/fsaverage_sym/label/s2fsaverage_sym_figures/Figure_W1"
      echo "Output (D): requires MRI to run, skipped."
      matlab -nodesktop -nosplash -r "figure_t14d1_raw;exit()"
      echo "Output (D): stats/figures/T14d1/Figure_W2-3A"
      matlab -nodesktop -nosplash -r "figure_t14d1;exit()"
      echo "Output (E): stats/figures/T14d1/Figure_W2-3B"
      break
      ;;
    "fig s9") #14 - running - done
      matlab -nodesktop -nosplash -r "figure_t11_new;exit()"
      echo "Output (A-C): stats/figures/T11_new"
      matlab -nodesktop -nosplash -r "s2fsaverage;exit()"
      echo "Output (D): data/coregistration/fsaverage_sym/label/s2fsaverage_sym_figures/Figure_W1"
      matlab -nodesktop -nosplash -r "figure_t14d1_raw;exit()"
      echo "Output (D): stats/figures/T14d1/Figure_W2-3A"
      matlab -nodesktop -nosplash -r "figure_t14d1;exit()"
      echo "Output (E): stats/figures/T14d1/Figure_W2-3B"
      break
      ;;
    "fig s10") #15
      matlab -nodesktop -nosplash -r "figure_t16_dk;exit()" #wait on brainexport script to finish
      echo "Output (A): stats/figures/T16"
      matlab -nodesktop -nosplash -r "figure_t19;exit()"
      echo "Output (B): stats/figures/T19"
      break
      ;;
    "fig s11") #16
      matlab -nodesktop -nosplash -r "figure_t14_allatl;exit()"
      echo "Output: stats/figures/T14_allatl"
      break
      ;;
    "fig s12") #17
      matlab -nodesktop -nosplash -r "fig_S16;exit()"
      echo "Output (A): stats/figures/figure_S16_pial_nbip-150_nedge-0_1stdev"
      matlab -nodesktop -nosplash -r "fig_cluster3_v2;exit()"
      echo "Output (B): stats/figures/cluster3_v2"
      matlab -nodesktop -nosplash -r "figure_t14_150;exit()"
      echo "Output (C): stats/figures/T14_150"
      break
      ;;
    "fig s13")
      matlab -nodesktop -nosplash -r "figure_t14_150;exit()"
      echo "Output (A): stats/figures/T14_150"
      matlab -nodesktop -nosplash -r "figure_t18;exit()"
      echo "Output (A): stats/figures/T18"
      break
      ;;
    "fig s14")
      matlab -nodesktop -nosplash -r "figure_t14_allatl;exit()"
      echo "Output (A-G): stats/figures/T14_allatl/atl2_Desikan-Killiany"
      matlab -nodesktop -nosplash -r "figure_t14_150;exit()"
      echo "Output (B-H): stats/figures/T14_allatl/atl2_Desikan-Killiany"
      break
      ;;
    quit)
      break
      ;;
    *) 
      echo "Invalid option $REPLY"
      ;;
  esac
done

