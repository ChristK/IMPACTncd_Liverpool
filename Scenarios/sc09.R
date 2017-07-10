#cmpfile("./Scenarios/sc09.R")
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
cat("Loading scenario 9...\n")
ptm <- proc.time()
#input
intervention_year <- 2011L
intervention_year2 <- 2017L
hc_annual_coverage <- 0.138 # from paula's report
hc_annual_uptake <- 0.323 # from paula's report
invitation_cost <- 5.11
participation_cost <- 13.28
kismet <- TRUE

# TODO atm cases prevented from HC do not considered eligible and are not expsed to intervention. This is very small underestimation of effectiveness

# Population wide effect 20% SSB tax from
# Briggs ADM, Mytton OT, Kehlbacher A, et al. Overall and income specific effect on 
# prevalence of overweight and obesity of 20% sugar sweetened drink tax in UK: econometric
# and comparative risk assessment modelling study. BMJ 2013;347:f6189. doi:10.1136/bmj.f6189
# table 8

intervention_year <- intervention_year - init.year
intervention_year2 <- intervention_year2 - init.year

# BMI 
POP[ssb.tax.effect.mc, on = c("age", "qimd"), bmival_sc09 := bmival + effect]
POP[!between(age, 16, 84) | year < intervention_year2, bmival_sc09 := bmival]
setkey(POP, id, year) # Do not remove
POP[, bmival.cvdlag := shift_byid(bmival_sc09, cvd.lag.mc, 20, id)] # lag exposure

# recalculate qrisk
set(POP, NULL, "qrisk2", 0)
POP[,
    qrisk2 := QRisk(
      age,
      sex,
      af,
      ra,
      kiddiag,
      bpmed,
      diab_type_I,# diab type 1
      diabtotr,
      bmival_sc09,
      origin,
      famcvd,
      tctohdl,
      omsysval,
      smoke_cat,
      townsend
    )]
POP[, qrisk_grp := cut(qrisk2, 
                       breaks = c(0, 0.1, 0.2, Inf), 
                       labels = c("qriskupto10", "qrisk10to20", "qrisk20+"), 
                       include.lowest = T, 
                       right = F)]

# select eligible population. If changes after participating to HC are permanent this have no bias
set(POP, NULL, "eligible_sc09", 0L)
POP[between(age, 40, 74) &
      chd.incidence == 0 &
      stroke.incidence == 0 &
      chd.mortality == 0 &
      stroke.mortality == 0 &
      other.mortality == 0 &
      undiag.diab == 0 &
      bpmed == 0 & # not known hypertensive
      af == 0 &
      ra == 0 &
      kiddiag == 0,
    eligible_sc09 := 1L]


# calculate invitation probability 
tt <- rep(1, yearstoproject - intervention_year)
x <- vector("list", 0)
for (i in 1:(yearstoproject - intervention_year))
{
  if (between(i, 1, 5))
  {
    if (i > 1) tt[i] <- (tt[i-1] - hc_annual_coverage)
    x$year[i] <- i + intervention_year - 1L
    x$prob[i] <- hc_annual_coverage / tt[i]
  } else {
    tt[i] <- tt[i-1]
    x$year[i] <- i + intervention_year - 1L
    x$prob[i] <- hc_annual_coverage / tt[i]
  }
}
setDT(x)

# calculate uptake probability 
uptake_weights <- fread("./Scenarios/HC_weights.csv")
uptake_weights <- melt(uptake_weights, 1:3, variable.name = "qrisk_grp", value.name = "weights")
tt <- uptake_weights[weights > 0, min(weights)]
uptake_weights[weights == 0, weights := tt]
uptake_weights[, `:=` (sex = factor_(sex),
                       qimd = ordered(qimd)
)]

# Estimate coverage
POP[x, `:=` (coverage = prob), on = "year"]
set(POP, NULL, "invited_sc09", 0L)
set(POP, NULL, "dice_coverage",  dice(nrow(POP)))
setkey(POP, id, year) # Do not remove
POP[eligible_sc09 == 1 & year >= intervention_year, 
    invited_sc09 := HC_coverage(id, coverage, dice_coverage)]

# Estimate uptake
set(POP, NULL, "participated_sc09", 0L)
set(POP, NULL, "dice_uptake",  dice(nrow(POP)))
POP[uptake_weights, `:=` (uptake = weights),
    on = c("agerange", "sex", "qimd", "qrisk_grp")]
POP[invited_sc09 == 1,
    uptake := hc_annual_uptake * uptake * .N,
    by = .(year)] #bayes theorem numerator
POP[invited_sc09 == 1,
    uptake := uptake / .N,
    by = .(agerange, sex, qimd, qrisk_grp, year)] #bayes theorem denominator
POP[invited_sc09 == 1, uptake := fcompress(uptake, .N * hc_annual_uptake), by=year]# Compress weights to givel reach uptake level
POP[invited_sc09 == 1, participated_sc09 := rbinom(.N, 1, uptake)]

# Estimate Prescription rate and persistence
prescription_weights <- fread("./Scenarios/HC_prescription.csv")
prescription_weights <- melt(prescription_weights, 1, variable.name = "qrisk_grp", value.name = "weights")
prescription_weights[, `:=` (qimd = ordered(qimd))]
POP[prescription_weights, prescription_uptake := weights, on = c("qimd", "qrisk_grp")]
POP[participated_sc09  == 1,  statins_uptake := bound(prescription_uptake * .N / count_if(cholval >= 5)), by = .(qimd, qrisk_grp)] # adjust weights to reflect only those with high cholesterol
POP[participated_sc09  == 1 & cholval >= 5,
    statins_uptake := fcompress(normalise(cholval, na.rm = T), sum(statins_uptake)),  
    by = .(qimd, qrisk_grp)] # give higher probabilities to those with higher chol  
POP[participated_sc09  == 1,  antihtn_uptake := bound(prescription_uptake * .N / count_if(omsysval >= 135)), by = .(qimd, qrisk_grp)] # adjust weights to reflect only those with high SBP
POP[participated_sc09  == 1 & omsysval >= 135,
    antihtn_uptake := fcompress(normalise(omsysval, na.rm = T), sum(antihtn_uptake)),  
    by = .(qimd, qrisk_grp)] # give higher probabilities to those with higher SBP  

set(POP, NULL, "statins_prescribed_sc09", 0L)
POP[participated_sc09  == 1 & cholval >= 5,  statins_prescribed_sc09 := rbinom(.N, 1, statins_uptake)]
set(POP, NULL, "statins_taken_sc09", 0L)
POP[statins_prescribed_sc09  == 1,  statins_taken_sc09 := rbinom(.N, 1, persistence.mc)]

set(POP, NULL, "antihtn_prescribed_sc09", 0L)
POP[participated_sc09  == 1 & omsysval >= 135,  antihtn_prescribed_sc09 := rbinom(.N, 1, antihtn_uptake)]
set(POP, NULL, "antihtn_taken_sc09", 0L)
POP[antihtn_prescribed_sc09  == 1,  antihtn_taken_sc09 := rbinom(.N, 1, persistence.mc)]

set(POP, NULL, "participation_effect_sc09", 0L)
setkey(POP, id, year) # Do not remove
POP[, participation_effect_sc09 := HC_effect(id, participated_sc09, invited_sc09)] # Effect lasts until next invitation
POP[, statins_taken_sc09 := med_effect(id, statins_taken_sc09, participation_effect_sc09)] # Effect lasts until next invitation
POP[, antihtn_taken_sc09 := med_effect(id, antihtn_taken_sc09, participation_effect_sc09)] # Effect lasts until next invitation
#View(POP[year >= intervention_year& eligible_sc09 == 1, .(id, year, invited_sc09, participated_sc09, participation_effect_sc09, statins_prescribed_sc09, statins_taken_sc09)])

# Estimate effect
#cholesterol
POP[, cholval_sc09 := cholval]
POP[statins_taken_sc09  == 1,  cholval_sc09 := cholval * (1 - atorv.eff.mc * adherence.mc)]
setkey(POP, id, year) # Do not remove
POP[, cholval.cvdlag := shift_byid(cholval_sc09, cvd.lag.mc, 3.8, id)] # lag exposure
#SBP
POP[, omsysval_sc09 := omsysval]
POP[antihtn_taken_sc09  == 1,  omsysval := omsysval - (omsysval - 135) * adherence.mc] # Assume that antihtn medication can potentialy achieve sbp 135 for all
POP[omsysval < omsysval, omsysval := omsysval]
setkey(POP, id, year) # Do not remove
POP[, omsysval.cvdlag := shift_byid(omsysval, cvd.lag.mc, 115, id)] # lag exposure

#recalculate diabetes prob
# Diab probabilities
POP[, bpmed := bpmed_init] #reset for each scenario. not needed for sc09. Here for consistency
POP[antihtn_taken_sc09 == 1L, bpmed := 1L] # feedback loop. Antihtn meds increased prescription increases diab 
set(POP, NULL, "pr_diab_nocvd", 0)
set(POP, NULL, "pr_diab_cvd", 0)
POP[between(age, 25, 84), c("pr_diab_nocvd", "pr_diab_cvd") := QDrisk(age, sex, b_corticosteroids, bpmed, bmival.cvdlag, origin, fh.diab, smoke_cat, townsend, 1L)]

#recalculate chd prob
set(POP, NULL, "chd.bmi.rr",  1)
POP[chd.bmi.rr.mc, chd.bmi.rr := bound(rr^((bmival.cvdlag - 20)/4.56), 1, 20), on = "age"]
set(POP, NULL, "chd.chol.rr",  1)
POP[chol.rr.chd.mc, chd.chol.rr := bound(rr^(3.8 - cholval.cvdlag), 1, 20), on = "age"]
set(POP, NULL, "chd.sbp.rr",  1)
POP[sbp.rr.chd.mc, chd.sbp.rr := bound(rr^((115 - omsysval.cvdlag)/20), 1, 20),
    on = c("age", "sex")]
set(POP, NULL, "stroke.bmi.rr",  1)
POP[stroke.bmi.rr.mc, stroke.bmi.rr := bound(rr^((bmival.cvdlag - 20) / 4.56), 1, 20),
    on = "age"]
set(POP, NULL, "stroke.chol.rr",  1)
POP[chol.rr.stroke.mc, stroke.chol.rr := bound(rr^(3.8 - cholval.cvdlag), 1, 20),
    on = "age"]
set(POP, NULL, "stroke.sbp.rr",  1)
POP[sbp.rr.stroke.mc, stroke.sbp.rr := bound(rr^((115 - omsysval.cvdlag)/20), 1, 20),
    on = c("age", "sex")]

set(POP, NULL, "pr_chd_inc", 0)
POP[between(age, ageL, ageH), pr_chd_inc := 
      p0_chd * chd.tob.rr * chd.ets.rr *
      chd.sbp.rr * chd.chol.rr * chd.pa.rr *
      chd.bmi.rr * chd.fv.rr]
set(POP, NULL, "pr_stroke_inc", 0)
POP[between(age, ageL, ageH), pr_stroke_inc := 
      p0_stroke * stroke.tob.rr * stroke.ets.rr *
      stroke.sbp.rr * stroke.chol.rr * stroke.pa.rr *
      stroke.bmi.rr * stroke.fv.rr]

set(POP, NULL, "htn.incidence_sc09", 0L)
POP[omsysval > 135 | bpmed == 1, htn.incidence_sc09 := 6000L]

if (kismet == F)
{
  POP[, `:=` (dice_chd_inc      = dice(.N), 
              dice_chd_death    = dice(.N),
              dice_stroke_inc   = dice(.N), 
              dice_stroke_death = dice(.N),
              dice_other_death  = dice(.N),
              dice_diab_inc     = dice(.N))] 
}

# Calculate incidence
setkey(POP, id, year)
POP[, diabtotr_sc09 := diabtotr_init]
set(POP, NULL, "chd.incidence_sc09", 0L)
set(POP, NULL, "stroke.incidence_sc09", 0L)
POP[year >= 0, c("diabtotr_sc09", "chd.incidence_sc09", "stroke.incidence_sc09") :=
      incidence_events(diabtotr_init, pr_diab_cvd, pr_diab_nocvd, dice_diab_inc,
                       chd.incidence_init, pr_chd_inc, chd.diab.rr, dice_chd_inc,
                       stroke.incidence_init, pr_stroke_inc, stroke.diab.rr,
                       dice_stroke_inc, id, year, cvd.lag.mc)]
POP[, diabtotr_sc09 := factor_(diabtotr_sc09)]

# Calculate mortality
set(POP, NULL, "death.tob.rr", 1)
POP[cigst1.cvdlag == "4", death.tob.rr := smokriskofdeath]
set(POP, NULL, "death.diab.rr", 1)
POP[, diabtotr.cvdlag := shift_byid(diabtotr_sc09, cvd.lag.mc, 1L, id, 1:2)] # lag exposure
POP[diabtotr.cvdlag == "2", death.diab.rr := 1.6] #rr from DECODE study

POP[, pr_other_death := qx]
POP[year >= 0 & between(age, ageL, ageH), pr_other_death :=
      bound(qx * death.diab.rr * death.tob.rr)] #

POP[, `:=` (chd.mortality_sc09 = 0L,
            stroke.mortality_sc09 = 0L,
            other.mortality_sc09 = 0L)]
setkey(POP, id, year)
POP[year >= 0 & chd.incidence_sc09 > 0,
    chd.mortality_sc09 := mortality_events(pr_chd_death, dice_chd_death, id, -20L)]
POP[year >= 0 & stroke.incidence_sc09 > 0,
    stroke.mortality_sc09 := mortality_events(pr_stroke_death, dice_stroke_death,
                                             id, -3L)]
POP[year >= 0,
    other.mortality_sc09 := mortality_events(pr_other_death, dice_other_death,
                                            id, -100L)]


# Output ------------------------------------------------------------------
out <- prepare_outcomes(POP)

# collect risk factors
RF <- collect_RF(POP)
nam <- paste0(init.year:(init.year + yearstoproject - 1))
fwrite_safe(summary_RF_fn(RF, nam, "sc09"),         output.dir("sc_RF.csv"),      append = TRUE)
fwrite_safe(summary_RF_fn(RF, nam, "sc09", "qimd"), output.dir("sc_RF_qimd.csv"), append = TRUE)
fwrite_safe(summary_RF_fn(RF, nam, "sc09", "sex"),  output.dir("sc_RF_sex.csv"),  append = TRUE)
fwrite_safe(summary_RF_fn(RF, nam, "sc09", "age"),  output.dir("sc_RF_age.csv"),  append = TRUE)

# shape output
out <- mask_age(out)

fwrite_safe(summary_fn(out, nam, "sc09"),         output.dir("sc_results.csv"),      append = TRUE)
fwrite_safe(summary_fn(out, nam, "sc09", "qimd"), output.dir("sc_results_qimd.csv"), append = TRUE)
fwrite_safe(summary_fn(out, nam, "sc09", "sex"),  output.dir("sc_results_sex.csv"),  append = TRUE)
fwrite_safe(summary_fn(out, nam, "sc09", "age"),  output.dir("sc_results_age.csv"),  append = TRUE)

POP[, grep("_sc09", names(POP), value = T) := NULL]
rm(out, nam, tt, x, RF)

gc()
print(proc.time() - ptm)
