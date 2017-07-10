#cmpfile("./life table engine.R")
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

# Load files --------------------------------------------------------------
cat("Generating life table...\n\n")
# Create diseasestoexclude.ICD from diseasestoexclude user input
diseasestoexclude.ICD <- diseasestoexclude
if ("CHD" %in% diseasestoexclude) {
  remove <- "CHD"
  diseasestoexclude.ICD <- outersect(diseasestoexclude.ICD, remove)
  diseasestoexclude.ICD <-
    c(diseasestoexclude.ICD, "Ischaemic heart diseases")
}

if ("stroke" %in% diseasestoexclude) {
  remove <- "stroke"
  diseasestoexclude.ICD <- outersect(diseasestoexclude.ICD, remove)
  diseasestoexclude.ICD <-
    c(diseasestoexclude.ICD,"Cerebrovascular diseases")
}

# if ("C16" %in% diseasestoexclude) {
#   remove <- "C16"
#   diseasestoexclude.ICD <- outersect(diseasestoexclude.ICD, remove)
#   diseasestoexclude.ICD <-
#     c(diseasestoexclude.ICD,"Malignant neoplasm of stomach")
# }
# 
# if ("C34" %in% diseasestoexclude) {
#   remove <- "C34"
#   diseasestoexclude.ICD <- outersect(diseasestoexclude.ICD, remove)
#   diseasestoexclude.ICD <-
#     c(diseasestoexclude.ICD,"Malignant neoplasm of trachea, bronchus and lung")
# }

population <-
  fread(
    "./Population/Population by adjusted IMD 2010 quintile_final.csv",
    header = T, 
    skip = 0
  )

deaths.qimd <-
  fread(
    "./LifeTables/Deaths by adjusted IMD 2010 quintile_final.csv",
    header = T, 
    skip = 0
  )

deaths.causes <-
  fread(
    "./LifeTables/deaths from selected causes by quintile_final_tcm77-388639.csv",
    header = T,
    skip = 0
  )

deaths.causes <- deaths.causes[cause %in% diseasestoexclude.ICD,]

deaths.causes <-
  deaths.causes[, list(disease = sum(disease)), 
                by = .(year, agegroup, sex, qimd)]

tt <-
  merge(
    population, deaths.qimd, by = c("year", "sex", "qimd", "agegroup"),
    all.x = T
  )
tt <-
  merge(
    tt, deaths.causes, by = c("year", "sex", "qimd", "agegroup"), 
    all.x = T
  )
tt[agegroup %!in% unique(agegroup.fn(ageL:ageH)), disease := 0] 
tt[agegroup %in% c("0", "01-04"), agegroup := "00-04"]
tt <- tt[, .(pop=sum(pop), deaths=sum(deaths), disease = sum(disease)),
         by = .(year, sex, qimd, agegroup)]
tt[, agegroup := relevel(factor(agegroup), "00-04")]
tt[, Mx.all := deaths / pop]  
#death rate all diseases tt[as.numeric(as.character(agegroup)) < ageL | 
#as.numeric(as.character(agegroup)) > ageH, disease := 0] # Exclude disease 
#mortality only for ages between ageL and ageH

# Exclude disease mortality only for ages between ageL and ageH
tt[, Mx.disease := (deaths - disease) / (pop - disease)] # from above only for ages 30-84
# death rate modelled diseases excluded

# file name for projection file
nam <- paste0("./LifeTables/", "mortality_projections_", paste(sort(diseasestoexclude), collapse = "_"), ".rds")

# Load projection from disk if exists ------------------------

if (file.exists(nam)) pp <- readRDS(nam) else {
  dependencies("demography")
  # Life tables loop --------------------------------------------------------
  
  # hor <- yearstoproject - 2013 + init.year
  # if (hor < 1) hor <- 1
  hor <- 60L # maximum simulation horizon
  
  # Full mortality projections (no disease excluded)
  if (Sys.info()[1] == "Linux")
    registerDoParallel(clusternumber)
  if (Sys.info()[1] == "Windows") {
    cl <- makeCluster(clusternumber)
    registerDoParallel(cl)
  }
  pp <- 
    foreach(
      k = 1:2,
      .inorder = F, 
      .combine = 'c'
    ) %:% # sex
    foreach(
      l = 1:5,
      .inorder = F,
      .packages = c("demography", "reshape2", "data.table")
    ) %dopar% {
      output <- vector("list", 3)
      # all cause mortality
      xx <- demogdata(
        acast(tt[sex == k & qimd == l,],
              agegroup ~ year, value.var = "Mx.disease"),
        acast(tt[sex == k & qimd == l,],
              agegroup ~ year, value.var = "pop"),
        c(seq(2.5, 87.5, 5), 92.5), #c(0:99),
        #c(0, 1, 5, seq(10, 85, 5), 100), 
        tt[sex == k & qimd == l, unique(year)],
        "mortality",
        paste0("England.sex-", k, ".qimd-", l), 0
      )
      xx <-
        smooth.demogdata(
          xx, age.grid = 0:100
        )
      temp <-
        forecast.fdm(fdm(xx,  method = "M", max.age = 100), hor)
      xx <- combine.demogdata(xx, temp)
      xx <- lifetable(xx, type = "period")
      xx <-
        as.data.table(xx$qx, keep.rownames = T)[, `:=`(sex = k, qimd = l)]
      setnames(xx, "rn", "age")
      xx[, `:=`(sex = factor(sex), age = as.numeric(age), qimd = ordered(qimd))]
      output[[1]] <- xx
      # minus disease specific mortality
      xx <- demogdata(
        acast(tt[sex == k & qimd == l,],
              agegroup ~ year, value.var = "Mx.all"),
        acast(tt[sex == k & qimd == l,],
              agegroup ~ year, value.var = "pop"),
        c(seq(2.5, 87.5, 5), 92.5), #c(0:99),
        #c(0, 1, 5, seq(10, 85, 5), 100),
        tt[sex == k &
             qimd == l, unique(year)],
        "mortality",
        paste0("England.sex-", k, ".qimd-", l), 0
      )
      #plot(xx)
      xx <-
        smooth.demogdata(
          xx, age.grid = 0:100
        )
      xxx <- xx$pop
      temp <-
        forecast.fdm(fdm(xx,  method = "M", max.age = 100), hor)
      
      xx <- combine.demogdata(xx, temp)
      #plot(xx)
      xx <- lifetable(xx, type = "period")
      
      xx <-
        as.data.table(xx$qx, keep.rownames = T)[, `:=`(sex = k, qimd = l)]
      xxx <-
        as.data.table(data.frame(xxx),
                      keep.rownames = T)[, `:=`(sex = k, qimd = l)]
      setnames(xx, "rn", "age")
      setnames(xxx, "rn", "age")
      setnames(xxx, paste0("X0.", 2002:2013), paste0(2002:2013))
      xx[, `:=`(sex = factor(sex), age = as.numeric(age), qimd = ordered(qimd))]
      xxx[, `:=`(sex = factor(sex), age = as.numeric(age), qimd = ordered(qimd))]
      output[[2]] <- xx
      output[[3]] <- xxx
      output
    }
  
  if (Sys.info()[1] == "Windows") stopCluster(cl)
  
  saveRDS(pp, nam, compress = FALSE)
  
  rm(hor)
}

# Lifetable.diab <- rbindlist(lapply(pp, `[[`, 2))
# Lifetable.diab[, age := as.integer(age)]
Lifetable <-
  rbindlist(lapply(pp, `[[`, 1)) # access first element of a nested list
Lifetable[, age := as.integer(age)]

# setkey(Lifetable.diab, age, sex, qimd)

Lifetable <- 
  melt(Lifetable, c("age", "sex", "qimd"), variable.name = "year", value.name = "qx", variable.factor = F)[, year := as.integer(year) - 2011L]
setkey(Lifetable, age, sex, qimd)

if (alignment) Lifetable[qimd == "5", qx := bound(1.17 * qx)]
if (alignment) Lifetable[qimd == "4", qx := bound(1.1 * qx)]

# Population distribution  ----------------------
# Match the sex and age structure of the initial year
population.actual <- fread("./Population/Population by adjusted IMD 2010 quintile_Liverpool.csv", header = T)
population.projections <- fread("./Population/Population projections_Liverpool.csv", header = T)

population.projections30 <- melt(
  population.projections[agegroup %in% unique(agegroup.fn(30:34))],
  1:2, variable.name = "year", value.name = "pop")[
    ,`:=`(pop = pop * 1e3, year = as.integer(factor_to_char(year)))]
population.projections30 <- rbind(
  population.actual[between(age, 30, 34), .(pop = sum(pop)), by = year],
  population.projections30[year > 2014, .(pop = sum(pop)), by = year]
)
population.projections30[, pop.growth := pop/population.projections30[year == init.year, pop]] 

population.projections <- melt(
  population.projections[agegroup %in% unique(agegroup.fn(30:84))],
  1:2, variable.name = "year", value.name = "pop")[
    ,`:=`(pop = pop * 1e3, year = as.integer(factor_to_char(year)))]
population.projections <- rbind(
  population.actual[between(age, 30, 84), .(pop = sum(pop)), by = year],
  population.projections[year > 2014, .(pop = sum(pop)), by = year]
)
#fwrite(population.projections, "./Population/population.projections30to84.csv")
population.projections[, pop.growth := pop/population.projections[year == init.year, pop]] 

# Calculate the exact fraction of the mid init.year population this sample represents
population.actual <- population.actual[year == init.year]
pop.fraction      <- n / population.actual[, sum(pop)]
sum_pop           <- population.actual[, sum(pop)]
                                       
cat(paste0("Population fraction = ", pop.fraction, "\n"), 
    file = output_dir("simulation parameters temp.txt"),
    append = T)

# Oversample population below 30 to approximate future population distribution
population.actual[between(age, 0, 30),
                  pop := population.actual[
                    sex == .BY[[1]] & qimd == .BY[[2]] & age == 30, pop],
                  by = .(sex, qimd)]
# population.actual[, year := bound(30 - age + init.year, 0, population.projections30[, max(year)])]
# population.actual[population.projections30[year > init.year],
#                   pop := as.integer(round(pop * pop.growth)),
#                   on = "year"]
# population.actual[, year := NULL]
population.actual[, pct := round(pop * n / sum_pop)]

# In Liverpool there is a huge peak of student population aged 19-24 that do not stay in the city after their studies.
# The following correction prevent them from being sampled initially, that would distort the future model population   
# population.actual[between(age, 19, 25), pct := population.actual[
#   (age == 18 | age == 26) & sex == .BY[[1]] & qimd == .BY[[2]], mean(pct)],
#   by = .(sex, qimd)]



# Smoking death-rate inflation tables -------------------------------------
# Load rr of death for smokers from Peto, Mortality in relation to smoking: 50
# years observations on male British doctors. Table 1
smokriskofdeath <- fread("./LifeTables/smokriskofdeath.csv")
setkey(smokriskofdeath)
smokriskofdeath <-
  smokriskofdeath[!diseasestoexclude,] 
# Need to adjust 1st column for each new disease
smokriskofdeath <-
  smokriskofdeath[, (sum(Current) * 4680 / 35.40) / (sum(Lifelong.non.smokers) *
                                                       2917 / 19.38)] 


rm(tt, deaths.qimd, population, deaths.causes, diseasestoexclude.ICD,
   pp, nam, population.projections30)
