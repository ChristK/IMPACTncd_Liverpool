#!/usr/bin/Rscript
## IMPACTncd: A decision support tool for primary prevention of NCDs Copyright (C)
## 2015 Chris Kypridemos

## IMPACTncd is free software; you can redistribute it and/or modify it under the
## terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.

## This program is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/> or write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301 USA.


# ********************************************************************************
# IMPACT NCD Liverpool
# ********************************************************************************
# User input
n <- 1e5  # Define the sample size
yearstoproject <- 30L  # NEED TO force >=1 and up to 50
numberofiterations <- 21L
clusternumber <- 20L # Change to your number of CPU cores
process_output <- T
clear_intermediate_outputs <- T # If T delete auxiliary output directories when simulation finish
init.year <- 2011L
ageL <- 30L  # Define lower age limit to diseases-model simulation (min = 30)
ageH <- 84L  # Define lower age limit to diseases-model simulation (max = 84)
widesynthpop <- F
alignment <- T # T or F (apply correction factor for lung cancer model)
qdrisk <- T # Use QDrisk score for diabetes incidence
Fertility.Assumption <- "N"  # Select (N)ormal, (H)igh or (L)ow fertility rate asumptions based on ONS scenarios. They do matter for accurate population statistics
cvd.lag <- 5L
fatality.annual.improvement.chd    <- 5L # 3 means 5% annual improvement in fatality
fatality.annual.improvement.stroke <- 2L # 3 means 2% annual improvement in fatality
diseasestoexclude <- c("CHD", "stroke")  # Define disease to be excluded from lifetables
# ICD10 code reminder for disease coding (http://apps.who.int/classifications/icd10/browse/2010/en#/I20-I25)

# *****************************************************************************
gc()
options(datatable.verbose = F)

cat("Initialising IMPACTncd...\n\n")
options(warn = 1)

get_dropbox_path <- function(x = character(0)) {
  # if (!require(jsonlite)) stop("You need to install jsonlite package.")
  if (!require(jsonlite)) install.packages("jsonlite")
  if (Sys.info()[1] == "Windows")
  {
    if (file.exists(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json")))
    {
      dropbox_path <- read_json(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json"))$personal$path
    } else dropbox_path <- read_json(paste0(Sys.getenv("LOCALAPPDATA"), "/Dropbox/info.json"))$personal$path
  } else {
    dropbox_path <- read_json("~/.dropbox/info.json")$personal$path
  }
  return(paste0(dropbox_path, "/", x))
}

if (grepl("incline", Sys.info()["nodename"], fixed = T))
{
  setwd("~/IMPACTncd Liverpool/")
  clusternumber <- ifelse(clusternumber > 30, 30, clusternumber)
  output_dir <- function(x = character(0))
    (paste0("./Output/", x))
  synthpop_dir <- "./SynthPop"
} else if (Sys.info()["nodename"] == "dh012142") {
  setwd(get_dropbox_path("Models/IMPACTncd Liverpool/"))
  clusternumber <- ifelse(clusternumber > 25, 25, clusternumber)
  output_dir <-
    function(x = character(0))
      (paste0("/mnt/storage_slow/Model_Outputs/Liverpool_Model/", x))
  synthpop_dir <-
    "~/Synthetic Populations/Northwest England 2011/"
} else {
  setwd(get_dropbox_path("Models/IMPACTncd Liverpool/"))
  output_dir <- function(x = character(0))
    (paste0("./Output/", x))
  synthpop_dir <- "./SynthPop"
}


source(file = "./initialisation.R")
source(file = "./life table engine.R")
source(file = "./diseases epidemiology.R")

# cl <- makeCluster(clusternumber) # used for clustering. win compatible
# registerDoParallel(cl)
registerDoParallel(clusternumber)  # used for forking. only linux
time.mark("start parallelisation")
foreach(iterations = 1:numberofiterations, .inorder = F, .verbose = T,
        .packages = c("data.table", "randtoolbox", "Kmisc", "compiler",
                      "mc2d", "Rcpp", "RcppArmadillo", "ckmisc"),
        .export = ls(), .noexport = c("scenarios.list", "time.mark")) %dorng%
  {
    my.env <- environment()  # get environment of this branch

    sys.source(file = "./cluster functions.R", my.env)
    sys.source(file = "./load synthetic population.R", my.env)

    sys.source(file = "./risk factor trajectories.R", my.env)
    sys.source(file = "./2dmc.R", my.env)

    time.mark("start simulation")
    sys.source(file = "./simulation.R", my.env)
    rm(my.env)  # BE CAREFULL. DTs altered with in my.env, change universaly.
    # You need copies of DTs to by handled within my.env
    return(0)
  }

if (exists("cl")) stopCluster(cl)
time.mark("End of parallelisation")
gc()

# Output
file.rename(output_dir("simulation parameters temp.txt"), output_dir("simulation parameters.txt"))
if (process_output == T) {
  source(file = "./post_simulation_functions.R")
  source(file = "./output.R")
}

end()

# compile scenarios
# lapply(list.files(path = './Scenarios', pattern = glob2rx('*.R'), full.names = T, recursive = F), cmpfile)
# lapply(list.files(path = './', pattern = glob2rx('*.R'), full.names = T, recursive = F), cmpfile)
