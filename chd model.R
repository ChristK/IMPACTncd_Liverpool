#cmpfile("./chd model.R")
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
cat("Loading CHD (I20-I25) model...\n")
ptm <- proc.time()
set(POP, NULL, "chd.incidence",  0L) 
set(POP, NULL, "chd.incidence_init",  0L) 

# RR for tobacco from Ezzati M, Henley SJ, Thun MJ, Lopez AD. Role of Smoking in Global and Regional 
# Cardiovascular Mortality. Circulation. 2005 Jul 26;112(4):489–97.
# Table 1 Model B

# RR for ex-smokers from Huxley RR, Woodward M. 
# Cigarette smoking as a risk factor for coronary heart disease
# in women compared with men: a systematic review and meta-analysis of prospective cohort studies. 
# The Lancet. 2011 Oct 14;378(9799):1297–305. 
# Appendix webfigure 8
#cat("smoking RR\n")
set(POP, NULL, "chd.tob.rr",  1)
POP[tobacco.rr.chd.mc, chd.tob.rr := rr, on = c("age", "sex", "cigst1.cvdlag")]

# RR for ETS He J, Vupputuri S, Allen K, Prerost MR, Hughes J, Whelton PK. Passive Smoking and the Risk of 
# Coronary Heart Disease — A Meta-Analysis of Epidemiologic Studies. New England Journal of Medicine. 1999;340(12):920–6. 
# Table 3. Adjusted RR
set(POP, NULL, "chd.ets.rr",  1)
POP[cigst1 %in% c("1", "2", "3") & expsmokCat.cvdlag == "1", chd.ets.rr := chd.ets.rr.mc]
POP[age > 69 & chd.ets.rr > 1, chd.ets.rr := bound(chd.ets.rr * (1-(age-69)/100), 1, 20)] # decrease risk for elderly

# RR for SBP from Optimal SBP level at 115mmHg and RR(HR) of dying from CHD was taken from "Age-specific relevance of
# usual blood pressure to vascular mortality: 
# a meta-analysis of individual data for one million adults in 61 prospective studies. 
# The Lancet. 2002 Dec 14;360(9349):1903–1913" 
# Figure 5
#cat("sbp RR\n")
set(POP, NULL, "chd.sbp.rr",  1)
POP[sbp.rr.chd.mc, chd.sbp.rr := bound(rr^((115 - omsysval.cvdlag)/20), 1, 20),
    on = c("age", "sex")]

# RR for Chol from "Blood cholesterol and 
# vascular mortality by age, sex, and blood pressure: a meta-analysis of individual data from 61 prospective studies 
# with 55.000 vascular deaths. The Lancet. 2007 Dec 7;370(9602):1829–39. 
# Appendix Webtable 6  fully adjusted
#cat("chol RR\n")
set(POP, NULL, "chd.chol.rr",  1)
POP[chol.rr.chd.mc, chd.chol.rr := bound(rr^(3.8 - cholval.cvdlag), 1, 20), on = "age"]

# RR for BMI from "The Emerging Risk Factors Collaboration.
# Separate and combined associations of body-mass index and abdominal adiposity with cardiovascular disease:
# collaborative analysis of 58 prospective studies.
# The Lancet 2011;377:1085–95. doi:10.1016/S0140-6736(11)60105-0
# Table 1 (Adjusted for age, sex, smoking status, systolic blood pressure, history of diabetes, and total and HDL cholesterol)
# and figure 2 for age specific gradient
#cat("bmi RR\n")
set(POP, NULL, "chd.bmi.rr",  1)
POP[chd.bmi.rr.mc, chd.bmi.rr := bound(rr^((bmival.cvdlag - 20)/4.56), 1, 20), on = "age"]

# RR for diab from The Emerging Risk Factors Collaboration. Diabetes mellitus, fasting blood glucose concentration, 
# and risk of vascular disease: a collaborative 
# meta-analysis of 102 prospective studies. The Lancet 2010;375:2215–22. doi:10.1016/S0140-6736(10)60484-9
# figure 2 (HRs were adjusted for age, smoking status, body-mass index, and  systolic blood pressure)
#cat("diab RR\n")
set(POP, NULL, "chd.diab.rr", 1)
POP[chd.diab.rr.mc, chd.diab.rr := rr, on = "age"] # The rr if diabetics. To be used later in random walk

# Specifically for PARF I need to calculate 
set(POP, NULL, "chd.diab.rr.forparf", 1)
POP[year == 0 & diabtotr == "2", chd.diab.rr.forparf := chd.diab.rr]

# RR for F&V from From Dauchet L, Amouyel P, Hercberg S, Dallongeville J. Fruit and Vegetable Consumption and Risk of Coronary Heart Disease: 
# A Meta-Analysis of Cohort Studies. J Nutr. 2006 Oct 1;136(10):2588–93. 
# To avoid negative PAF an optimal level of F&V has to be set arbitrarily. I set it to 10 
# when convert porftvg from categorical to numeric I create bias. eg 1=less than 1 portion
#cat("fv RR\n")
set(POP, NULL, "chd.fv.rr",  1)
POP[porftvg.cvdlag < 8L, chd.fv.rr := bound(chd.fv.rr.mc^(porftvg.cvdlag - 7L), 1, 20)] # max effect at 7 portions
POP[age > 69 & chd.fv.rr > 1, chd.fv.rr := bound(chd.fv.rr * (1-(age-69)/100), 1, 20)] # decrease risk for elderly

# RR for PA 1. WHO | Comparative Quantification of Health Risks [Internet]. 
# WHO [cited 2014 Jan 30];Available from: http://www.who.int/publications/cra/en/
# Table 10.19 (with adjustment for measurement error)
#cat("pa RR\n")
set(POP, NULL, "chd.pa.rr",  1)
POP[pa.rr.chd.mc, chd.pa.rr := rr, on = c("age", "a30to06m.cvdlag")]

# Estimate prevalence -----------------------------------------------------
#cat(paste0("Estimating CHD prevalence in ", init.year, " ...\n\n"))
age.structure <- setkey(POP[age <= ageH & year == 0, .N, by = .(age, sex, qimd)], age, sex, qimd)
age.structure[chd_epi$prevalence[between(age, 1, ageH)],
              Nprev := rbinom(.N, ifelse(is.na(N), 0L, N), prevalence), 
              on = c("age", "sex", "qimd")]
age.structure[chd_epi$incidence[between(age, 1, ageH)], 
              Nprev := ifelse(is.na(Nprev), 0L, Nprev) -
                rbinom(.N, ifelse(is.na(N), 0L, N) - ifelse(is.na(Nprev),
                                                            0L, Nprev), incidence),
              on = c("age", "sex", "qimd")]
age.structure[Nprev < 0, Nprev := 0]
setnames(age.structure, "N", "population")


id.chd <- POP[between(age, 1, ageH) & year == 0, 
              .(id = resample(id, age.structure[age == .BY[[1]] & sex == .BY[[2]] & qimd == .BY[[3]],
                                       Nprev], 
                     replace = F,
                     prob = chd.tob.rr * chd.ets.rr * 
                       chd.sbp.rr * chd.chol.rr * chd.bmi.rr * 
                       chd.diab.rr.forparf * chd.fv.rr * chd.pa.rr 
                     )), 
              by = .(age, sex, qimd, year)]
id.chd[chd_epi$duration, duration := as.integer(round(bound(rnorm(.N, chd.duration, se), 1, Inf))), on = "age"] # estimate duration
id.chd[age - duration < 0, duration := age]
POP[id.chd, chd.incidence_init := duration, on = c("id", "year")]

rm(id.chd)

# Estimate PAF ------------------------------------------------------------
#cat("Estimating CHD PAF...\n")
  chdpaf <- 
    POP[between(age, ageL, ageH) & chd.incidence_init == 0L & year == 0, 
        .(paf = 1 - 1 / (sum(chd.tob.rr * chd.ets.rr *
                               chd.sbp.rr * chd.chol.rr * 
                               chd.bmi.rr * chd.diab.rr.forparf * chd.fv.rr *
                               chd.pa.rr) / .N)), 
        by = .(age, sex, qimd)
        ]
  setkey(chdpaf, age, sex, qimd)
  chdpaf[, paf := predict(loess(paf~age, span = 0.5)), by = .(sex, qimd)]
  # chdpaf[, {
  #   plot(age, paf, main = paste0(.BY[[1]],"-", .BY[[2]]))
  #   lines(age, paf2)
  # }
  # , keyby = .(sex, qimd)]
  setkey(chd_epi$incidence, age, sex, qimd)
  chd_epi$incidence[chdpaf, p0 := incidence * (1 - paf)]
  chd_epi$incidence[is.na(p0), p0 := incidence]


POP[chd_epi$incidence, p0_chd := p0, on = c("age", "sex", "qimd")]

# Estimate CHD incidence -------------------------------
#cat("Estimating CHD incidence without diabetes...\n\n")
set(POP, NULL, "pr_chd_inc", 0)
POP[between(age, ageL, ageH), pr_chd_inc := 
      p0_chd * chd.tob.rr * chd.ets.rr *
      chd.sbp.rr * chd.chol.rr * chd.pa.rr *
      chd.bmi.rr * chd.fv.rr]

POP[chd_epi$fatality, pr_chd_death := fatality * ((100-fatality.annual.improvement.chd)/100)^year,
    on = c("age", "sex", "qimd")]
print(proc.time() - ptm)



