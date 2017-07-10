#cmpfile("./stroke model.R")
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
cat("Loading stroke (I60-I69) model...\n")
ptm <- proc.time()
set(POP, NULL, "stroke.incidence",  0L) 
set(POP, NULL, "stroke.incidence_init",  0L) 

# RR for tobacco from Ezzati M, Henley SJ, Thun MJ, Lopez AD. Role of Smoking in Global and Regional 
# Cardiovascular Mortality. Circulation. 2005 Jul 26;112(4):489–97.
# Table 1 Model B
#cat("smoking RR\n")
set(POP, NULL, "stroke.tob.rr",  1)
POP[tobacco.rr.stroke.mc, stroke.tob.rr := rr, on = c("age", "sex", "cigst1.cvdlag")]

#ex-smokers
# Cigarette smoking as a risk factor for stroke: The Framingham study
# "Stroke risk decreased significantly by two years and was at the
# level of nonsmokers by five years after cessation of cigarette smoking"

# Calculate PAF of ETS for stroke
# RR from Oono IP, Mackay DF, Pell JP. Meta-analysis of the association between secondhand smoke exposure and stroke. 
# J Public Health 2011;33:496–502. doi:10.1093/pubmed/fdr025
set(POP, NULL, "stroke.ets.rr",  1)
POP[cigst1 %in% c("1", "2", "3") & expsmokCat.cvdlag == "1", stroke.ets.rr := stroke.ets.rr.mc]
POP[age > 69 & stroke.ets.rr > 1, stroke.ets.rr := bound(stroke.ets.rr * (1-(age-69)/100), 1, 20)] # decrease risk for elderly

# Calculate RR for stroke. Optimal SBP level at 115mmHg and RR(HR) of dying from
# stroke was taken from "Age-specific relevance of usual blood pressure to 
# vascular mortality: a meta-analysis of individual data for one million adults in 61 prospective studies. 
# The Lancet. 2002 Dec 14;360(9349):1903–1913" 
# Figure 3
#cat("sbp RR\n")
set(POP, NULL, "stroke.sbp.rr",  1)
POP[sbp.rr.stroke.mc, stroke.sbp.rr := bound(rr^((115 - omsysval.cvdlag)/20), 1, 20),
    on = c("age", "sex")]

# Calculate RR for stroke. Optimal chol level at 3.8 mmol/L and RR(HR) of 
# dying from stroke was taken from "Blood cholesterol and 
# vascular mortality by age, sex, and blood pressure: a meta-analysis of 
# individual data from 61 prospective studies 
# with 55.000 vascular deaths. The Lancet. 2007;370:1829–39. 
# Figure 4 (for total stroke. I used only significant HR's).
#cat("chol RR\n")
set(POP, NULL, "stroke.chol.rr",  1)
POP[chol.rr.stroke.mc, stroke.chol.rr := bound(rr^(3.8 - cholval.cvdlag), 1, 20),
    on = "age"]

# RR for BMI from "The Emerging Risk Factors Collaboration.
# Separate and combined associations of body-mass index and abdominal adiposity
# with cardiovascular disease: collaborative analysis of 58 prospective studies.
# The Lancet 2011;377:1085–95. doi:10.1016/S0140-6736(11)60105-0
# Table 1 (Adjusted for age, sex, smoking status, systolic blood pressure,
# history of diabetes, and total and HDL cholesterol)
# BMI not significant for ischaemic stroke but other obesity metrics are. 
#!! NEED TO decide if I want to use it
#cat("bmi RR\n")
set(POP, NULL, "stroke.bmi.rr",  1)
POP[stroke.bmi.rr.mc, stroke.bmi.rr := bound(rr^((bmival.cvdlag - 20) / 4.56), 1, 20),
    on = "age"]

# RR for diabetes from The Emerging Risk Factors Collaboration. 
# Diabetes mellitus, fasting blood glucose concentration, 
# and risk of vascular disease: a collaborative 
# meta-analysis of 102 prospective studies. The Lancet 2010;375:2215–22
# figure 2 (HRs were adjusted for age, smoking status, body-mass index, 
# and  systolic blood pressure)
#cat("diab RR\n")
set(POP, NULL, "stroke.diab.rr",  1)
POP[stroke.diab.rr.mc, stroke.diab.rr := rr, on = "age"]# The rr if diabetics. To be used later in random walk

# Specifically for PARF I need to calculate 
set(POP, NULL, "stroke.diab.rr.forparf",  1)
POP[year == 0 & diabtotr == "2", stroke.diab.rr.forparf := stroke.diab.rr]

# Calculate RR for stroke. From Dauchet L, Amouyel P, Dallongeville J. 
# Fruit and vegetable consumption and risk of stroke A meta-analysis of cohort 
# studies. Neurology. 2005;65:1193–7. 
# To avoid negative PAF an optimal level of F&V has to be set arbitrarily. I set it to 8 
#cat("fv RR\n")
set(POP, NULL, "stroke.fv.rr",  1)
POP[porftvg.cvdlag < 8L, stroke.fv.rr := stroke.fv.rr.mc^(porftvg.cvdlag - 7L)] 
POP[age > 69 & stroke.fv.rr > 1, stroke.fv.rr := bound(stroke.fv.rr * (1-(age-69)/100), 1, 20)] # decrease risk for elderly

# RR for PA 1. WHO | Comparative Quantification of Health Risks [Internet]. 
# WHO [cited 2014 Jan 30];Available from: http://www.who.int/publications/cra/en/
# Table 10.20 (with adjustment for measurement error)
#cat("pa RR\n")
set(POP, NULL, "stroke.pa.rr",  1)
POP[pa.rr.stroke.mc, stroke.pa.rr := rr, on = c("age", "a30to06m.cvdlag")]

# Estimate prevalence -----------------------------------------------------
age.structure <- setkey(POP[age <= ageH & year == 0, .N, by = .(age, sex, qimd)], age, sex, qimd)
age.structure[stroke_epi$prevalence[between(age, 1, ageH)],
              Nprev := rbinom(.N, ifelse(is.na(N), 0L, N), prevalence), 
              on = c("age", "sex", "qimd")]
age.structure[stroke_epi$incidence[between(age, 1, ageH)], 
              Nprev := ifelse(is.na(Nprev), 0L, Nprev) -
                rbinom(.N, ifelse(is.na(N), 0L, N) - ifelse(is.na(Nprev),
                                                            0L, Nprev), incidence),
              on = c("age", "sex", "qimd")]
age.structure[Nprev < 0, Nprev := 0]
setnames(age.structure, "N", "population")


id.stroke <- POP[between(age, 1, ageH) & year == 0, 
              .(id = resample(id, age.structure[age == .BY[[1]] & sex == .BY[[2]] & qimd == .BY[[3]],
                                                Nprev], 
                              replace = F,
                              prob = stroke.tob.rr * stroke.ets.rr * 
                                stroke.sbp.rr * stroke.chol.rr * stroke.bmi.rr * 
                                stroke.diab.rr.forparf * stroke.fv.rr * stroke.pa.rr 
              )), 
              by = .(age, sex, qimd, year)]
id.stroke[stroke_epi$duration, duration := as.integer(round(bound(rnorm(.N, stroke.duration, se), 1, Inf))), on = "qimd"] # estimate duration
id.stroke[age - duration < 0, duration := age]
POP[id.stroke, stroke.incidence_init := duration, on = c("id", "year")]

rm(id.stroke)

# Estimate PAF ------------------------------------------------------------
#cat("Estimating stroke PAF...\n")
strokepaf <- 
  POP[between(age, ageL, ageH) & stroke.incidence_init == 0L & year == 0, 
      .(paf = 1 - 1 / (sum(stroke.tob.rr * stroke.ets.rr *
                             stroke.sbp.rr * stroke.chol.rr * 
                             stroke.bmi.rr * stroke.diab.rr.forparf * stroke.fv.rr *
                             stroke.pa.rr) / .N)), 
      by = .(age, sex, qimd)
      ]
setkey(strokepaf, age, sex, qimd)
strokepaf[, paf := predict(loess(paf~age, span = 0.5)), by = .(sex, qimd)]
# strokepaf[, {
#   plot(age, paf, main = paste0(.BY[[1]],"-", .BY[[2]]))
#   lines(age, paf2)
# }
# , keyby = .(sex, qimd)]
setkey(stroke_epi$incidence, age, sex, qimd)
stroke_epi$incidence[strokepaf, p0 := incidence * (1 - paf)]
stroke_epi$incidence[is.na(p0), p0 := incidence]


POP[stroke_epi$incidence, p0_stroke := p0, on = c("age", "sex", "qimd")]

# Estimate stroke incidence -------------------------------
#cat("Estimating stroke incidence without diabetes...\n\n")
set(POP, NULL, "pr_stroke_inc", 0)
POP[between(age, ageL, ageH), pr_stroke_inc := 
      p0_stroke * stroke.tob.rr * stroke.ets.rr *
      stroke.sbp.rr * stroke.chol.rr * stroke.pa.rr *
      stroke.bmi.rr * stroke.fv.rr]

POP[stroke_epi$fatality, pr_stroke_death := fatality * ((100-fatality.annual.improvement.stroke)/100)^year,
    on = c("age", "sex", "qimd")]
print(proc.time() - ptm)
