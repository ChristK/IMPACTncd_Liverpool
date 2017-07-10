#cmpfile("./load synthetic population.R")
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


# Define function for output dir

cat("Loading synthetic population...\n\n")

SPOP <- readRDS(file = random.spop.file[[iterations]])

setDT(SPOP)


  SPOP[, `:=`(hsize             = NULL, 
              diabtyper         = NULL, 
              bpmedc            = NULL,
              lipid             = NULL,
              segment           = NULL,
              saltCat.intersalt = NULL,
              smokyrs           = as.integer(smokyrs)
              )]  


SPOP[cigst1 == "1", smokyrs := 0L]
# These were not able to quit for at least a year so then the relapse would not apply for them
SPOP[cigst1 %in% c("2","3") & endsmoke == 0, 
     `:=` (cigst1 = "4", cigdyal = numsmok, numsmok = 0)] 
SPOP[cigst1 == "4" & cigdyal < 1, 
     `:=` (cigdyal = 1)] 

agegroup.fn(SPOP)

# Stratified sampling
SPOP[age>90, age := 90]
SPOP[, id:= seq_len(.N)] # id is also .I in SPOP
tt <- SPOP[age > 0, sample(id, population.actual[age==.BY[[1]] & sex==.BY[[2]] & qimd==.BY[[3]], pct], T), by = .(age, sex, qimd)][, V1]
POP <- SPOP[tt, ] # because id is also .I in SPOP

POP[, id:= sample(.N)] # randomly allocate ids.

# datasets for ageing.distr function
bmi.rank <- setkey(
  SPOP[bmival>0, .(bmival, 
                      "bmival.rank" = perc.rank(bmival)), by = group],
  group, bmival.rank)

sbp.rank <- setkey(
  SPOP[omsysval>0, .(omsysval, 
                        "omsysval.rank" = perc.rank(omsysval)), by = group],
  group, omsysval.rank)

chol.rank <- setkey(
  SPOP[cholval>0, .(cholval, 
                       "cholval.rank" = perc.rank(cholval)), by = group],
  group, cholval.rank)


# Define origin (ethnicity) -----------------------------------------------
tt <- fread("./Population/Population by adjusted IMD 2010 quintile_Liverpool_ethnicity.csv", 
            colClasses = c("character", "integer", "numeric"))
tt[, qimd := ordered(qimd)]
tt[, race := as.integer(substring(race, 1, 1))]
tt[, cum.pr := cumsum(pr), by = qimd]
tt <- dcast(tt, qimd~race, value.var = "cum.pr")

# 1 = white # 2 = indian # 3 = pakistani # 4 = bangladeshi # 5 = other asian
# 6 = black caribbean # 7 = black african # 8 = chinese # 9 = other
tt <- POP[, .(id, qimd)][tt, on = "qimd"]
set(tt, NULL, "d", dice(nrow(tt)))
for (k in 3:11) set(tt, NULL, k, tt[, Reduce(`<`, .SD), .SDcol = c(k, 12)])
tt[,origin :=  1L + Reduce(`+`, .SD), .SDcol = 3:11]
POP[tt, origin := origin,
on = "id"]


# Estimate Townsend score
POP[qimd == "1", townsend := runif(.N,   -7, -3.4)]
POP[qimd == "2", townsend := runif(.N, -3.4, 0.2)]
POP[qimd == "3", townsend := runif(.N,  0.2, 3.8)]
POP[qimd == "4", townsend := runif(.N,  3.8, 7.4)]
POP[qimd == "5", townsend := runif(.N,  7.4, 11)]

POP[cigst1 == "4" & cigdyal < 1L, cigdyal := 1L]
POP[cigst1 == "3" & numsmok < 1L, numsmok := 1L]

POP[, `:=` (cigdyal = as.numeric(cigdyal), numsmok = as.numeric(numsmok))]

index <- POP[, .SD, .SDcols = c("id", "age", "sex", "qimd")]
setkey(index, id)

rm(tt, SPOP)
# sink() 
