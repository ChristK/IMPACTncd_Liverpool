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

# preample ----------------------------------------------------------------
gc()
dependencies <- function(x) {
  for (j in x) {
    # require returns T invisibly if it was able to load package
    if (!require(j, character.only = T)) {
      # If package was not able to be loaded then re-install
      install.packages(j, dependencies = T)
      # Load package after installing
      require(j, character.only = T)
    }
  }
}



# Then try/install packages...
dependencies(c("simPopulation",
               "data.table",
               "Hmisc",
               "doParallel",
               "compiler",
               "survey",
               "MASS",
               "glm2",
               "StatMatch",
               "pryr",
               "mc2d",
               "rriskDistributions",
               "quantreg",
               "foreach",
               "rms"))



# Get Dropbox folder.  Automaticly define of my dropbox folder in windows Manually define on linux
options(warn = 1)

if (Sys.info()[1] == "Linux") {
  if (system("whoami", T )== "mdxasck2") {
    setwd("~/IMPACTncd Liverpool/")
    clusternumber <- ifelse (clusternumber < 30, 30, clusternumber)  # overwrites previous if <60
  } else {
    setwd(paste("/home/",
                system("whoami", T),
                "/Dropbox/Models/IMPACTncd Liverpool/",
                sep = "",
                collapse = ""))
    get.dropbox.folder <- function() paste0("/home/", system("whoami", T),
                                            "/Dropbox/")
  }
} else if (Sys.info()[1] == "Darwin") {
  setwd("/Users/chris/Dropbox/Models/IMPACTncd Liverpool/")
} else {
  get.dropbox.folder <- function() {
    if (!require(RCurl))
      stop("You need to install RCurl package.")
    if (Sys.info()["sysname"] != "Windows")
      stop("Currently, 'get.dropbox.folder' works for Windows and Linux only. Sorry.")
    db.file <- paste(Sys.getenv("LOCALAPPDATA"), "\\Dropbox\\host.db", sep = "")
    base64coded <- readLines(db.file, warn = F)[2]
    base64(base64coded, encode = F)
  }
  setwd(paste0(get.dropbox.folder(), "/Models/IMPACTncd Liverpool/"))
}



# Define function to clear labels form SPSS labelled imports
clear.labels <- function(x) {
  if (is.list(x)) {
    for (i in 1 : length(x)) class(x[[i]]) <- setdiff(class(x[[i]]), "labelled")
    for (i in 1 : length(x)) attr(x[[i]],"label") <- NULL
  }
  else {
    class(x) <- setdiff(class(x), "labelled")
    attr(x, "label") <- NULL
  }
  return(invisible(x))
}

# Define function to calculate SD from svy objects (not necessary)
# svysd <- function(...) sqrt(coef(svyvar(...)))

# from plyr. relevel vector
mapvalues <- function (x, from, to, warn_missing = TRUE)
{
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }
  if (is.factor(x)) {
    levels(x) <- mapvalues(levels(x), from, to, warn_missing)
    return(x)
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ",
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
  }
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}

# Define function to split agegroups and create groups
agegroup.fn <- function(x, lagtime = 0) {
  breaks                   <- c(0, 1, seq(5, 85, 5), 130)
  labels                   <- c("<1   ", "01-04", "05-09",
                                "10-14", "15-19", "20-24",
                                "25-29", "30-34", "35-39",
                                "40-44", "45-49", "50-54",
                                "55-59", "60-64", "65-69",
                                "70-74", "75-79", "80-84", "85+")
  if (is.numeric(x)) {
    agegroup = cut(x - lagtime,
                   breaks = breaks,
                   labels = labels,
                   include.lowest = T,
                   right = F,
                   ordered_result = T)
    return(invisible(agegroup))
  } else {
    if (is.data.table(x)) {
      x[, agegroup := cut(as.numeric(as.character(age)) + lagtime,
                          breaks = breaks,
                          labels = labels,
                          include.lowest = T,
                          right = F,
                          ordered_result = T)]
      x[, group := paste(qimd, sex, agegroup, sep = '')]
      return(invisible(x))
    } else return(print("not eligible input"))
  }
}
options(survey.lonely.psu = "adjust") #Lonely PSU (center any single-PSU strata around the sample grand mean rather than the stratum mean)
#options(datatable.auto.index= F) # needed until bug is resolved

# ************************************************************** CHD PREVALENCE *****************************************************************
load(file = paste0(get.dropbox.folder(), "/Datasets/Health Survey for England/2011/hse2011ai.RData"))
HSE2011original <- clear.labels(HSE.2011)
rm(HSE.2011)
HSE2011original <- data.table(HSE2011original, key="age")
agegroup.fn(HSE2011original)
#HSE0.srv.blood <- svydesign(id=~psu, strata=~cluster, weights=~wt.blood, nest=F, data=HSE2011original[wt.blood>0], check.strata = T)
#HSE0.srv.nurse <- svydesign(id=~psu, strata=~cluster, weights=~wt.nurse, nest=F, data=HSE2011original[wt.nurse>0], check.strata = T)
#HSE0.srv.int <- svydesign(id=~psu, strata=~cluster, weights=~wt.int, nest=F, data=HSE2011original[wt.nurse>0], check.strata = T)

HSE2011 = copy(HSE2011original[, list(wt.int, wt.nurse, wt.blood, psu, cluster, age, agegroup, sex, group, qimd, bmival, cholval1, omsysval, diabtotr, cigst1, startsmk, endsmoke, porftvg, frtpor, numsmok, smokyrs, cigdyal, origin, hdlval1,
                                      iregdef, bpmedd, diabete2, sha, ihdis2, strodef, ageheart, ageangi, agestro)])
HSE2011[, `:=`(year=0, a30to06 = NA, sodium = NA, potass = NA, creatin = NA, wt.urine = 0)]
HSE2011[, psu := paste0(psu, "2011")]
HSE2011[, cluster := paste0(cluster, "2011")]
HSE2011[, sha := as.integer(sha)]
agegroup.fn(HSE2011)


HSE2011[, table(ihdis2)]
HSE2011[, table(strodef)]

# I will use IHDIS2: (D) Had IHD (Angina or Heart Attack) {revised} 1 = Yes, 2  = No (Apparently the same as ihdis)
HSE2011[is.na(ihdis2) == T, ihdis2 := 99] # Convert NAs to 99 = NA
HSE2011[, ihdis2 := factor(ihdis2)]
HSE2011[is.na(strodef) == T, strodef := 99] # Convert NAs to 99 = NA
HSE2011[, strodef := factor(strodef)]
HSE.srv <- svydesign(id=~psu, strata=~cluster, weights=~wt.int, nest=F, data=HSE2011)

CHD <- data.table(svyby(~ihdis2, ~agegroup + sex + qimd, HSE.srv, svymean))
CHD[, `:=` (ihdis22 = NULL, se.ihdis22 = NULL, ihdis299 = NULL, se.ihdis299 = NULL, se.ihdis21 = NULL)]
setnames(CHD, "ihdis21", "prevalence")
CHD[, prevalence.sm := predict(loess(prevalence~as.integer(agegroup))), by = .(sex, qimd)]
CHD[prevalence.sm < 0, prevalence.sm := 0 ]
CHD[, fwrite(.SD, paste0("./CVD Statistics/dismod/chd prevalence_sex_", .BY[[1]], "_qimd_", .BY[[2]], ".csv")), by = .(sex, qimd)]

stroke <- data.table(svyby(~strodef, ~agegroup + sex + qimd, HSE.srv, svymean))
stroke[, `:=` (strodef2 = NULL, se.strodef2 = NULL, strodef99 = NULL, se.strodef99 = NULL, se.strodef1 = NULL)]
setnames(stroke, "strodef1", "prevalence")
stroke[, prevalence.sm := predict(loess(prevalence~as.integer(agegroup))), by = .(sex, qimd)]
stroke[prevalence.sm < 0, prevalence.sm := 0 ]

stroke[, fwrite(.SD, paste0("./CVD Statistics/dismod/stroke prevalence_sex_", .BY[[1]], "_qimd_", .BY[[2]], ".csv")), by = .(sex, qimd)]

population <-
  fread(
    "./Population/Population by adjusted IMD 2010 quintile_final.csv",
    header = T,
    skip = 0
  )
population[year == 2011, fwrite(.SD, paste0("./CVD Statistics/dismod/pop_sex_", .BY[[1]], "_qimd_", .BY[[2]], ".csv")), by = .(sex, qimd)]


deaths.qimd <-
  fread(
    "./LifeTables/Deaths by adjusted IMD 2010 quintile_final.csv",
    header = T,
    skip = 0
  )
deaths.qimd[population, pop := pop, on = c("year", "agegroup", "sex", "qimd")]
deaths.qimd[, rate := deaths/pop]
deaths.qimd[year == 2011, fwrite(.SD, paste0("./CVD Statistics/dismod/alldeaths_sex_", .BY[[1]], "_qimd_", .BY[[2]], ".csv")), by = .(sex, qimd)]

deaths.causes <-
  fread(
    "./LifeTables/deaths from selected causes by quintile_final_tcm77-388639.csv",
    header = T,
    skip = 0
  )
deaths.causes[population, pop := pop, on = c("year", "agegroup", "sex", "qimd")]
deaths.causes[, rate := disease/pop]
deaths.causes[year == 2011 & cause == "Ischaemic heart diseases", fwrite(.SD, paste0("./CVD Statistics/dismod/chddeaths_sex_", .BY[[1]], "_qimd_", .BY[[2]], ".csv")), by = .(sex, qimd)]
deaths.causes[year == 2011 & cause == "Cerebrovascular diseases", fwrite(.SD, paste0("./CVD Statistics/dismod/strokedeaths_sex_", .BY[[1]], "_qimd_", .BY[[2]], ".csv")), by = .(sex, qimd)]

# Duration of diseases ------
HSE2011[ihdis2 == 1, chd.duration := age - ageangi]
HSE2011[ihdis2 == 1 & is.na(chd.duration), chd.duration := age - ageheart]
HSE2011[strodef == 1, stroke.duration := age - agestro]

HSE.srv <- svydesign(id=~psu, strata=~cluster, weights=~wt.int, nest=F, data=HSE2011)
HSE.srv <- subset(HSE.srv, ihdis2 == 1)

mm <- svyglm(chd.duration ~ agegroup + sex + qimd, HSE.srv )
anova(mm) # only age important
xx <- data.table(svyby(~chd.duration, ~ agegroup, HSE.srv, svymean, na.rm = T))

xx[, fwrite(.SD, paste0("./CVD Statistics/chd duration.csv"))]

HSE.srv <- svydesign(id=~psu, strata=~cluster, weights=~wt.int, nest=F, data=HSE2011)
HSE.srv <- subset(HSE.srv, strodef == 1)
mm <- svyglm(stroke.duration ~ qimd, HSE.srv )
anova(mm) #  qimd important, age but not agegroup also important
xx <- data.table(svyby(~stroke.duration, ~ qimd, HSE.srv, svymean, na.rm = T))

xx[, fwrite(.SD, paste0("./CVD Statistics/stroke duration.csv"))]
