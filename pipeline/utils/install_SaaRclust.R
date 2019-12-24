#!/usr/bin/env Rscript

#This Rscript installs SaaRclust package from the github repository along with all dependencies.
#author: David Porubsky

Sys.setenv(Renv='PWD')
library(devtools)

withr::with_libpaths(new = "utils/R-packages", install_git("git://github.com/maryamghr/SaaRclust.git", branch = "master"), "prefix")
