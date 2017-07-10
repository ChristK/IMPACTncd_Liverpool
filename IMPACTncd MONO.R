## IMPACTncd: A decision support tool for primary prevention of NCDs
## Copyright (C) 2015  Chris Kypridemos

## IMPACTncd is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>
## or write to the Free Software Foundation, Inc., 51 Franklin Street,
## Fifth Floor, Boston, MA 02110-1301  USA.


# Preample ----------------------------------------------------------------
n <- 5e4  # Define the sample size
yearstoproject <- 10L  # NEED TO force >=1 and up to 50
numberofiterations <- 5L
clusternumber <- 4L # Change to your number of CPU cores 
init.year <- 2011L
ageL <- 30L  # Define lower age limit to diseases-model simulation (min = 30)
ageH <- 84L  # Define lower age limit to diseases-model simulation (max = 84)
widesynthpop <- F
alignment <- T # T or F (apply correction factor for lung cancer model)
qdrisk <- T # Use QDrisk score for diabetes incidence
Fertility.Assumption <- "N"  # Select (N)ormal, (H)igh or (L)ow fertility rate asumptions based on ONS scenarios. They do matter for accurate population statistics
cvd.lag <- 5L 
fatality.annual.improvement.chd    <- 5L # 3 means 3% annual improvement in fatality
fatality.annual.improvement.stroke <- 2L # 2 means 2% annual improvement in fatality
diseasestoexclude <- c("CHD", "stroke")  # Define disease to be excluded from lifetables
# ICD10 code reminder for disease coding (http://apps.who.int/classifications/icd10/browse/2010/en#/I20-I25)
clear_intermediate_outputs <- F # If T delete auxiliary output directories when simulation finish
process_output <- T
# *****************************************************************************

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

# Main --------------------------------------------------------------------
gc()
options(datatable.verbose = F)

source(file = "./initialisation.R")
source(file = "./life table engine.R")
source(file = "./diseases epidemiology.R")

iterations = 2L 
my.env <- environment() # get environment of this branch
# Define functions in foreach loop
sys.source(file = "./cluster functions.R", my.env)
#sys.source(file = "./diseases epidemiology.R", my.env)
# Load synthetic population
sys.source(file = "./load synthetic population.R", my.env)

# Actual simulation
i = init.year - 2011L
sys.source(file = "./risk factor trajectories.R", my.env)
sys.source(file = "./2dmc.R", my.env)

sys.source(file = "./lifecourse engine.R", my.env)
sys.source(file = "./chd model.R", my.env)
sys.source(file = "./stroke model.R", my.env)

POP[, `:=` (dice_chd_inc = dice(.N), 
            dice_chd_death = dice(.N),
            dice_stroke_inc = dice(.N), 
            dice_stroke_death = dice(.N),
            dice_other_death = dice(.N),
            dice_diab_inc = dice(.N))] 

# Calculate incidence
setkey(POP, id, year)
POP[year >= 0, c("diabtotr", "chd.incidence", "stroke.incidence") :=
      incidence_events(diabtotr_init, pr_diab_cvd, pr_diab_nocvd, dice_diab_inc,
                       chd.incidence_init, pr_chd_inc, chd.diab.rr, dice_chd_inc,
                       stroke.incidence_init, pr_stroke_inc, stroke.diab.rr,
                       dice_stroke_inc, id, year, cvd.lag.mc)]

# Calculate mortality
sys.source(file = "./other model.R", my.env)
POP[, `:=` (chd.mortality = 0L,
            stroke.mortality = 0L,
            other.mortality = 0L)]
setkey(POP, id, year)
POP[year >= 0 & chd.incidence > 0,
    chd.mortality := mortality_events(pr_chd_death, dice_chd_death, id, -20L)]
POP[year >= 0 & stroke.incidence > 0,
    stroke.mortality := mortality_events(pr_stroke_death, dice_stroke_death,
                                         id, -3L)]
POP[year >= 0,
    other.mortality := mortality_events(pr_other_death, dice_other_death,
                                        id, -100L)]

sys.source(file = "./partial_output.R", my.env)

# Scenarios
sys.source(file = "./Scenarios/HC_scenario_prep.R", my.env)
lapply(scenarios.list, sys.source, envir = my.env)

if (process_output == T) {
  source(file = "./post_simulation_functions.R")
  source(file = "./output.R")
} 

# compile scenarios
#lapply(list.files(path = "./Scenarios", pattern = glob2rx("*.R"), full.names = T, recursive = F), cmpfile)
#lapply(list.files(path = "./", pattern = glob2rx("*.R"), full.names = T, recursive = F), cmpfile)
