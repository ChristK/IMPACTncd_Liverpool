#cmpfile("./diseases epidemiology.R")
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

# CHD statistics ----------------------------------------------------------
if ("CHD" %in% diseasestoexclude) {
  xx <-  rbind(
    fread("./CVD Statistics/CHD DISMOD Females 2011 QIMD1.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 1)],
    fread("./CVD Statistics/CHD DISMOD Females 2011 QIMD2.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 2)],
    fread("./CVD Statistics/CHD DISMOD Females 2011 QIMD3.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 3)],
    fread("./CVD Statistics/CHD DISMOD Females 2011 QIMD4.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 4)],
    fread("./CVD Statistics/CHD DISMOD Females 2011 QIMD5.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 5)],
    fread("./CVD Statistics/CHD DISMOD Males 2011 QIMD1.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 1)],
    fread("./CVD Statistics/CHD DISMOD Males 2011 QIMD2.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 2)],
    fread("./CVD Statistics/CHD DISMOD Males 2011 QIMD3.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 3)],
    fread("./CVD Statistics/CHD DISMOD Males 2011 QIMD4.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 4)],
    fread("./CVD Statistics/CHD DISMOD Males 2011 QIMD5.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 5)]    
  )[, `:=` (sex = factor(sex), age = as.integer(Age), qimd = as.ordered(qimd))]
  
  chd_epi <- vector("list", 0)
  chd_epi$incidence <- setnames(copy(xx[, c(15, 13, 14, 5), with = F]),
                                "Incidence (rates)", "incidence")[
                                  , incidence := as.numeric(incidence)]
  
  chd_epi$prevalence <- setnames(copy(xx[, c(15, 13, 14, 6), with = F]),
                                 "Prevalence (rates)", "prevalence")[
                                   , prevalence := as.numeric(prevalence)]
  
  chd_epi$fatality <- setnames(copy(xx[, c(15, 13, 14, 8), with = F]),
                               "Case fatality (rates)", "fatality")[
                                 , fatality := as.numeric(fatality)]
  if (alignment) chd_epi$fatality[qimd == "5", fatality :=
                                      bound(fatality * 1.05)]
  
  xx <- fread("./CVD Statistics/chd duration.csv")
  xxx <- data.table(age = 1:ageH)[, agegroup := agegroup.part(age)]
  xx <- xx[xxx, on = "agegroup"]
  xx[, agegroup := agegroup.fn(age)]
  xx[age < 40, `:=` (chd.duration = xx[age == 40, chd.duration], se = xx[age == 40, se])]
  rm(xxx)
  
  chd_epi$duration <- copy(xx)
  }


# Stroke statistics -------------------------------------------------------
# Do I have to separate between ischaemic and haemorrhagic? The risk factors seems more ore less the same.
if ("stroke" %in% diseasestoexclude) {
  xx <-  rbind(
    fread("./CVD Statistics/Stroke DISMOD Females 2011 QIMD1.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 1)],
    fread("./CVD Statistics/Stroke DISMOD Females 2011 QIMD2.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 2)],
    fread("./CVD Statistics/Stroke DISMOD Females 2011 QIMD3.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 3)],
    fread("./CVD Statistics/Stroke DISMOD Females 2011 QIMD4.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 4)],
    fread("./CVD Statistics/Stroke DISMOD Females 2011 QIMD5.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 2, qimd = 5)],
    fread("./CVD Statistics/Stroke DISMOD Males 2011 QIMD1.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 1)],
    fread("./CVD Statistics/Stroke DISMOD Males 2011 QIMD2.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 2)],
    fread("./CVD Statistics/Stroke DISMOD Males 2011 QIMD3.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 3)],
    fread("./CVD Statistics/Stroke DISMOD Males 2011 QIMD4.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 4)],
    fread("./CVD Statistics/Stroke DISMOD Males 2011 QIMD5.csv",
          sep = ",", header = T, stringsAsFactors = F,
          skip = 3, nrows = 85)[, `:=`(sex = 1, qimd = 5)]    
  )[, `:=` (sex = factor(sex), age = as.integer(Age), qimd = as.ordered(qimd))]
  
  stroke_epi <- vector("list", 0)
  stroke_epi$incidence <- setnames(copy(xx[, c(15, 13, 14, 5), with = F]),
                                "Incidence (rates)", "incidence")[
                                  , incidence := as.numeric(incidence)]
  
  stroke_epi$prevalence <- setnames(copy(xx[, c(15, 13, 14, 6), with = F]),
                                 "Prevalence (rates)", "prevalence")[
                                   , prevalence := as.numeric(prevalence)]
  
  stroke_epi$fatality <- setnames(copy(xx[, c(15, 13, 14, 8), with = F]),
                               "Case fatality (rates)", "fatality")[
                                 , fatality := as.numeric(fatality)]
  if (alignment) stroke_epi$fatality[qimd == "5", fatality :=
                                      bound(fatality * 1.12)]
  
  stroke_epi$duration <- fread("./CVD Statistics/stroke duration.csv")[
    , qimd := ordered(qimd)]
}


if (exists("xx")) rm(xx)

