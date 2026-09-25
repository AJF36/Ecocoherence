##################################################
# misosoup_analysis.R
##################################################
# In this script I perform ...
#
# Madrid, August 2026
# Spanish National Center for Biotechnology
# https://github.com/AJF36/
##################################################

rm(list=ls())

results_dir = "results/misosoup_analysis.R"
figures_dir = "figures/misosoup_analysis.R"

results_dir = c(results_dir,figures_dir)

for (dir in results_dir) {
if(!(fs::dir_exists(dir))) {,
fs::dir_create(dir) }
}
