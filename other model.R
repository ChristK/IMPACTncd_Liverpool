#cmpfile("./other model.R")
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
cat("Estimating deaths from other causes...\n")
ptm <- proc.time()
#cat("Inflate mortality for diabetics and smokers...\n\n")
#Doll R, et al. Mortality in relation to smoking: 50 years’ observations on male
#British doctors. BMJ 2004;328:1519. doi:10.1136/bmj.38142.554479.AE table 1
set(POP, NULL, "death.tob.rr", 1)
POP[cigst1.cvdlag == "4", death.tob.rr := smokriskofdeath]

set(POP, NULL, "death.diab.rr", 1)
POP[, diabtotr.cvdlag := shift_byid(diabtotr, cvd.lag.mc, 1L, id, 1:2)] # lag exposure
POP[diabtotr.cvdlag == "2", death.diab.rr := 1.6] #rr from DECODE study

deathpaf <- 
    POP[between(age, ageL, ageH) & year >= 0, #For consistency with scenarios. Otherwise it double counts smoking decrease and prevalence increase that already taken into account during forecasting
        .(paf = 1 - 1 / (sum(death.tob.rr * death.diab.rr) / .N)), 
        keyby = .(age, sex, qimd, year)
        ]
deathpaf[, paf := predict(loess(paf~age, span = 0.50)), keyby = .(sex, qimd, year)]

setkey(deathpaf, age, sex, qimd, year)
Lifetable[, qx_adj := qx]
Lifetable[deathpaf, qx_adj := qx * (1 - paf), on = c("age", "sex", "qimd", "year")]

POP[Lifetable, qx := qx_adj, on = c("age", "sex", "qimd", "year")]
POP[age > 100, qx := 1]
POP[year < 0, qx := 0]
POP[, pr_other_death := qx]
POP[year >= 0 & between(age, ageL, ageH), pr_other_death :=
      bound(qx * death.diab.rr * death.tob.rr)] #
print(proc.time() - ptm)
