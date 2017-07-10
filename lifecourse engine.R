#cmpfile("./lifecourseengine.R")
## IMPACTncd: A decision support tool for primary prevention of NCDs
## Copyright (C) 2016  Chris Kypridemos

## IMPACTncdUS is free software; you can redistribute it and/or modify
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
ptm <- proc.time()
setkey(POP, group)
POP[, a30to06m.rank := perc.rank(a30to06m),
    by = .(group)
    ] # this can overwrite a30to06m for efficiency
POP[is.na(a30to06m.rank), a30to06m.rank := dice(.N)] # fix for small groups
POP[, porftvg.rank := perc.rank(porftvg),
    by = .(group)
    ]
POP[is.na(porftvg.rank), porftvg.rank := dice(.N)]
POP[, bmival.rank := perc.rank(bmival),
    by = .(group)
    ]
POP[is.na(bmival.rank), bmival.rank := dice(.N)]
POP[, omsysval.rank := perc.rank(omsysval),
    by = .(group)
    ]
POP[is.na(omsysval.rank), omsysval.rank := dice(.N)]
POP[, cholval.rank := perc.rank(cholval),
    by = .(group)
    ]
POP[is.na(cholval.rank), cholval.rank := dice(.N)]
POP[cigst1 == "4", cigdyal.rank := perc.rank(cigdyal),
    by = .(group)
    ]
POP[is.na(cigdyal.rank), cigdyal.rank := dice(.N)] # for those who wll become smokers
POP[, diabtotr_init := diabtotr]

# On corticosteroids -------------
# randomly assign .5% of the population
POP[, `:=` 
    (b_corticosteroids = rbinom(.N, 1, 0.005))]

# Family history of diabetes ------------
# If the prob of being diab is ~8%. So the probability of having at least
# one of 3 family members with diabetes is 1 - (1-0.08)^3.
# I let family members vary between 2 and 2+rpois(n, 1)
# This does not account for the future increase of diabetes prevalence. 
# Therefore it underestimates. I multiply current diab prev with 1.2 to 
# adjust for that

POP[, fh.diab := 
      rbinom(.N, 1, 1 - (1 - (1.2 * sum(diabtotr == "2") / .N))^(2 + rpois(.N, 1)))]
print("before clone DT")
print(proc.time() - ptm)
ptm <- proc.time()
# clone DT as many times as years + max(lag)
POP <- rbindlist(
  rep(
    list(POP),
    yearstoproject + cvd.lag.mc),#, cancer.lag)),
  idcol = "year"
)
POP[, year := year - cvd.lag.mc - 1L] # To be updated for cancers to max(cvd.lag.mc, cancer.lag)
# advance age
POP[, age := age + year] 
# birth engine needs to go here

POP <- POP[age >= 0] # negative ages have been created also
agegroup.fn(POP) # this need to ran after frank so ranks
#remain stable through the life course 
setkey(POP, age, sex, qimd)
print("Clone DT")
print(proc.time() - ptm)
ptm <- proc.time()

# PA estimation
POP[between(age, 20, 84),
    a30to06m := pred.pa(year, age, sex, qimd, a30to06m.rank)]

POP[, a30to06m.rank := NULL]
#POP[between(age, 20, 84)][id == sample(id, 1), plot(age, a30to06m, ylim=c(0,7))]
#POP[between(age, 20, ageH), mean(a30to06m), keyby = year][, plot(year, V1, ylim=c(0,7))] # PA mean falling rapidly
print("PA estimation")
print(proc.time() - ptm)
ptm <- proc.time()
# FV estimation
POP[between(age, ageL, ageH),
    porftvg := pred.fv(year, age, sex, qimd, porftvg.rank)]
POP[, porftvg.rank := NULL]
#POP[between(age, ageL, ageH)][id == sample(id, 1), plot(age, porftvg, ylim=c(0,8))]
#POP[between(age, ageL, ageH), mean(porftvg), keyby = year][, plot(year, V1, ylim=c(0,8))] # FV mean falling rapidly
print("FV estimation")
print(proc.time() - ptm)
ptm <- proc.time()

# Smoking back projection
POP[cigst1 %in% c("2", "3") & year < 0, endsmoke := endsmoke + year]
POP[cigst1 %in% c("2", "3") & year < 0 & endsmoke <= 0, 
    `:=` (cigst1 = "4", endsmoke = 0, cigdyal = numsmok, numsmok = 0)]

POP[cigst1 == "4" & year < 0, smokyrs := smokyrs + year]
POP[cigst1 == "4" & year < 0 & smokyrs <= 0, 
    `:=` (cigst1 = "1", endsmoke = 0, cigdyal = 0, numsmok = 0, smokyrs = 0)]
print("Smoking back projection")
print(proc.time() - ptm)
ptm <- proc.time()
# Forward projection
# baseline prevalence at age 20
POP[age == 18 & year > 0, tmp := pred.sm0prev(year, age, sex, qimd)]
POP[tmp == 1, `:=` (cigst1 = "4", 
                    endsmoke = 0,
                    cigdyal = pred.cigdyal(year, age, sex, qimd, 1, dice(.N)), 
                    numsmok = 0, 
                    smokyrs = rpois(.N, 3))]

POP[age == 18 & year > 0 & cigst1 == "1", tmp := pred.exsm0prev(year, age, sex, qimd)]
POP[tmp == 1, `:=` (cigst1 = "3", 
                    endsmoke = 1L + rpois(.N, 1),
                    cigdyal = 0, 
                    numsmok = pred.cigdyal(year, age, sex, qimd, 1, dice(.N)), 
                    smokyrs = rpois(.N, 2))]
POP[, tmp := NULL]

# Generate smoking related probabilities
pr_relapse <- expand.grid(endsmoke = 1:10, sex = factor(1:2), qimd = ordered(1:5))
pr_relapse <- cbind.dt(pr_relapse, 1 - pred.ex0sm1(pr_relapse$endsmoke, pr_relapse$sex, pr_relapse$qimd, "log"))
pr_relapse <- df2mat(dcast(pr_relapse, sex+qimd~endsmoke, value.var = "V1"))

POP[age > 18 & year >= 0, `:=` (
  pr_init = 0,
  pr_cess = pred.sm0ex1(year, age, sex, qimd),
  dice = dice(.N)
)]

POP[between(age, 19, 20) & year >= 0, 
    `:=` (pr_init = pred.nev0sm1(year, age, sex))]

POP[cigst1 == "2", cigst1 := "3"]

# smoking random walk
setkey(POP, id, year) # sort by id
tt <- setDT(POP[age > 18L & year >= 0L, simsmok(cigst1, pr_init, pr_cess, 
                                                dice, id, age, sex, qimd,
                                                endsmoke, smokyrs, cigdyal,
                                                numsmok, year, cigdyal.rank, 
                                                pr_relapse, 6L)]) # last fn input the cutoff for relapse. ie no relapse after 6 yrs of successful cessation

tt[cigst1 == "4", cigdyal := pred.cigdyal(year, age, sex, qimd, smokyrs, cigdyal.rank)]
tt[, numsmok := numsmok_fix(cigst1, id, cigdyal, numsmok)]
# tt[id == sample(id, 1),{
#   par(mfrow=c(4,1))
#   plot(year, cigst1, ylim=c(1,4), main = unique(id))
#   plot(year, cigdyal, ylim=c(1,30), main = unique(id))
#   plot(year, smokyrs, ylim=c(0,40), main = unique(id))
#   plot(year, endsmoke, ylim=c(0,40), main = unique(id))
#   }]
POP[tt, on = c("id", "year"), 
    `:=` (cigst1 = i.cigst1,
          endsmoke = i.endsmoke,
          smokyrs = i.smokyrs,
          cigdyal = i.cigdyal,
          numsmok = i.numsmok)]
POP[, `:=` (pr_init = NULL,
            pr_cess = NULL,
            dice = NULL,
            cigdyal.rank = NULL)]

# needed for QRisk and QDrisk
POP[, smoke_cat := 0L]
POP[cigst1 == "3", smoke_cat := 1L]
POP[cigst1 == "4", smoke_cat := 3L]
POP[cigst1 == "4" & cigdyal < 10, smoke_cat := 2L]
POP[cigst1 == "4" & cigdyal > 19, smoke_cat := 4L]

rm(tt, pr_relapse)
#POP[between(age, ageL, ageH) & year >= 0][id == sample(id, 1), plot(age, as.integer(cigst1), ylim=c(0,4))]
# POP[between(age, ageL, ageH), sum(cigst1 == "1")/.N, keyby = .(year)][, plot(year, V1, ylim=c(0,1), main = "Never")]
# POP[between(age, ageL, ageH), sum(cigst1 == "3")/.N, keyby = year][, plot(year, V1, ylim=c(0,1), main = "Ex smokers")]
# POP[between(age, ageL, ageH), sum(cigst1 == "4")/.N, keyby = year][, plot(year, V1, ylim=c(0,1), main = "Current smokers")]
print("Smoking forward projection")
print(proc.time() - ptm)
ptm <- proc.time()

# ETS ---------------------------------------------------------------------
# assumes linear relation with 
# smoking prevalence by qimd, stratified by age and sex
pred.ets <-
  POP[year == 0, sm.pr := prop_if(cigst1 =="4"), by = .(qimd)][
    year == 0, ets.pr := prop_if(expsmokCat == "1"), by = .(agegroup, sex, qimd)][
      year == 0, .(ets.coef = mean(ets.pr/sm.pr)), by = .(agegroup, sex, qimd)]
POP[, sm.pr := sum(cigst1 =="4") / .N, by = .(year, qimd)]
POP[pred.ets[agegroup %in% unique(agegroup.fn(ageL:ageH))], 
    expsmokCat := as.character(rbinom(.N, 1, sm.pr * ets.coef)),
    on = c("agegroup", "sex", "qimd")]
POP[, c("ets.pr", "sm.pr") := NULL]
rm(pred.ets)
# POP[between(age, ageL, ageH), sum(expsmokCat == "1") / .N, keyby=year][, plot(year, V1)]
print("ETS")
print(proc.time() - ptm)
ptm <- proc.time()

# BMI estimation
ageing.distr(bmi.rank, "bmival")
POP[between(age, 20, ageH),
    bmival := bmival - mean(bmival),
    by = .(group, year)]
POP[between(age, 20, ageH),
    bmival := bmival +
      pred.bmi(year, age, sex, qimd, a30to06m)]
# POP[id==47654, plot(year, bmival, ylim=c(10,50), main=min(age))]
# POP[id==sample(id, 1), plot(age, bmival, ylim=c(10,50), main=min(id))]
# POP[between(age, ageL, ageH), mean(bmival), keyby=year][, plot(year, V1)]
print("BMI")
print(proc.time() - ptm)
ptm <- proc.time()

#SBP
ageing.distr(sbp.rank, "omsysval")
POP[between(age, 20, ageH),
    omsysval := omsysval - mean(omsysval),
    by = .(group, year)]
POP[between(age, 20, ageH),
    omsysval := omsysval +
      pred.sbp(year, age, sex, qimd, bmival, cigst1, a30to06m)]
# POP[id==47654, plot(year, omsysval, ylim=c(70,200), main=min(age))]
# POP[id==sample(id, 1), plot(age, omsysval, ylim=c(70,200), main=min(id))]
# POP[between(age, ageL, ageH), mean(omsysval), keyby=year][, plot(year, V1)]
print("SBP")
print(proc.time() - ptm)
ptm <- proc.time()

#TC
ageing.distr(chol.rank, "cholval")
POP[between(age, ageL, ageH),
    cholval := cholval - mean(cholval),
    by = .(group, year)]
POP[between(age, ageL, ageH),
    cholval := cholval +
      pred.chol(year, age, sex, qimd, bmival)]
# POP[id==47654, plot(year, cholval, ylim=c(2,12), main=min(age))]
# POP[id==sample(id, 1), plot(age, cholval, ylim=c(2,12), main=min(id))]
# POP[between(age, ageL, ageH), mean(cholval), keyby=year][, plot(year, V1)]
POP[, c("bmival.rank", "omsysval.rank", "cholval.rank") := NULL]
print("TC")
print(proc.time() - ptm)
ptm <- proc.time()

# BP meds prevalence 
set(POP, NULL, "bpmed", 0L)
POP[between(age, 20, ageH), bpmed := pred.bpmed(.N, year, age, sex, qimd, omsysval)] # 25 for QDrisk
POP[, bpmed_init := bpmed] #so each scenario can alter bpmed based on bpmed_init and 
# QDrisk keep working always on bpmed
# POP[id==sample(id, 1), plot(age, bpmed, ylim=c(-0.1,1.1), main=min(id))]
# POP[between(age, ageL, ageH), sum(bpmed==1, na.rm = T)/.N, by = age][, plot(age, V1)]
# POP[between(age, ageL, ageH), sum(bpmed==1, na.rm = T)/.N, by = round(omsysval)][, plot(round, V1)]
set(POP, NULL, "htn.incidence", 0L)
POP[omsysval > 135 | bpmed == 1, htn.incidence := 6000L]

print("BPmeds")
print(proc.time() - ptm)
ptm <- proc.time()

# TC to HDL estimation ----------------------------------------------------------
set(POP, NULL, "tctohdl", 0)
POP[between(age, ageL, ageH),
    tctohdl := pred.tctohdl(cholval, age, sex, qimd, bmival, a30to06m, cigst1)]
# POP[id==sample(id, 1), plot(age, tctohdl, ylim=c(0,10), main=min(id))]
# POP[between(age, ageL, ageH), mean(tctohdl, na.rm = T), by = age][, plot(age, V1)]
# POP[between(age, ageL, ageH), mean(tctohdl, na.rm = T), by = round(cholval)][, plot(round, V1)]
print("TCtoHDL")
print(proc.time() - ptm)
ptm <- proc.time()

# FamCVD ------------------------------------------------------------------
set(POP, NULL, "famcvd", 0L)
POP[between(age, ageL, ageH), famcvd := pred.famcvd(.N, age, qimd)]
print("FamCVD")
print(proc.time() - ptm)
ptm <- proc.time()

# AF prevalence ------------------------------------------------------------------
set(POP, NULL, "af", 0L)
POP[between(age, ageL, ageH), af := pred.af(.N, age, omsysval)]
print("AF")
print(proc.time() - ptm)
ptm <- proc.time()

# CKD 4/5 prevalence -------------------------------------------------------------
set(POP, NULL, "kiddiag", 0L)
POP[between(age, ageL, ageH), kiddiag := pred.kidfailgp(.N, age, omsysval)]
print("ckd")
print(proc.time() - ptm)
ptm <- proc.time()

# Rheum arthr incidence --------------------------------------------------
set(POP, NULL, "ra", 0L)
POP[RApreval.rr.l, on = c("age", "sex"), ra := rbinom(.N, 1, rr)]
print("Reum")
print(proc.time() - ptm)
ptm <- proc.time()

# Create lag time variables and smoothing

setkey(POP, id, year) # Do not remove
POP[, cigst1.cvdlag := shift_byid(cigst1, cvd.lag.mc, 1L, id, 1:4)] # lag exposure
# POP[, cigst1.cvdlag2 := shift(cigst1, cvd.lag.mc, 1L), by =  id]
# identical(POP$cigst1.cvdlag1, POP$cigst1.cvdlag2)
POP[, expsmokCat.cvdlag := shift_byid(expsmokCat, cvd.lag.mc, 1L, id)] # lag exposure
# POP[, expsmokCat.cvdlag2 := shift(expsmokCat, cvd.lag.mc, 1L), by =  id]
# identical(POP$expsmokCat.cvdlag1, POP$expsmokCat.cvdlag2)
POP[, omsysval2 := roll_mean(omsysval, 4), by = id] # smoothing
POP[, omsysval := omsysval2 * mean(omsysval) / mean(omsysval2), by = .(group, year)] # ensure group means remain the same after smoothing
POP[, omsysval2 := NULL]
POP[, omsysval.cvdlag := shift_byid(omsysval, cvd.lag.mc, 115, id)] # lag exposure
# POP[, omsysval.cvdlag2 := shift(omsysval, cvd.lag.mc, 115), by =  id]
# identical(POP$omsysval.cvdlag1, POP$omsysval.cvdlag2)
POP[, cholval2 := roll_mean(cholval, 4), by = id] # smoothing
POP[, cholval := cholval2 * mean(cholval) / mean(cholval2), by = .(group, year)] # ensure group means remain the same after smoothing
POP[, cholval2 := NULL]
POP[, cholval.cvdlag := shift_byid(cholval, cvd.lag.mc, 3.8, id)] # lag exposure
# POP[, cholval.cvdlag2 := shift(cholval, cvd.lag.mc, 3.8), by =  id]
# identical(POP$cholval.cvdlag1, POP$cholval.cvdlag2)
POP[, bmival2 := roll_mean(bmival, 4), by = id] # smoothing
POP[, bmival := bmival2 * mean(bmival) / mean(bmival2), by = .(group, year)] # ensure group means remain the same after smoothing
POP[, bmival2 := NULL]
POP[, bmival.cvdlag := shift_byid(bmival, cvd.lag.mc, 20, id)] # lag exposure
# POP[, bmival.cvdlag2 := shift(bmival, cvd.lag.mc, 20), by =  id]
# identical(POP$bmival.cvdlag1, POP$bmival.cvdlag2)
POP[, porftvg.cvdlag := shift_byid(porftvg, cvd.lag.mc, 0, id)] # lag exposure
# POP[, porftvg.cvdlag2 := shift(porftvg, cvd.lag.mc, 0), by =  id]
# identical(POP$porftvg.cvdlag1, POP$porftvg.cvdlag2)
POP[, a30to06m.cvdlag := shift_byid(a30to06m, cvd.lag.mc, 0, id)] # lag exposure
# POP[, a30to06m.cvdlag2 := shift(a30to06m, cvd.lag.mc, 0), by =  id]
# identical(POP$a30to06m.cvdlag1, POP$a30to06m.cvdlag2)
print("Lagtimes")
print(proc.time() - ptm)
ptm <- proc.time()

# Diab probabilities
set(POP, NULL, "pr_diab_nocvd", 0)
set(POP, NULL, "pr_diab_cvd", 0)
POP[between(age, 25, 84), c("pr_diab_nocvd", "pr_diab_cvd") := QDrisk(age, sex, b_corticosteroids, bpmed, bmival.cvdlag, origin, fh.diab, smoke_cat, townsend, 1L)]
print("Diabetes")
print(proc.time() - ptm)
ptm <- proc.time()

# setkey(POP, id, year)
# POP[, bmival2 := roll_mean(bmival, 4), by = id]
# POP[id==sample(id, 1),
#     {
#       plot(age, bmival, ylim=c(10,40), main=min(id))
#       lines(age, bmival2, ylim=c(10,40), col="red")
# } 
# ]
# POP[between(age, 20, 84), .(mean(bmival) / mean(bmival2)), by = year]
# POP[, bmival3 := bmival2 * mean(bmival) / mean(bmival2), by = .(group, year)]
# POP[between(age, 20, 84), .(mean(bmival) / mean(bmival3)), by = year]
# POP[id==sample(id, 1),
#     {
#       plot(age, bmival, ylim=c(10,40), main=min(id))
#       lines(age, bmival2, ylim=c(10,40), col="red")
#       lines(age, bmival3, ylim=c(10,40), col="blue")
#       
#     } 
#     ]

