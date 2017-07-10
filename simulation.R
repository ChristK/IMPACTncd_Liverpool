#cmpfile("./simulation.R")
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


# Year 0 is mid-year 2011 (or any other init.year) Dates of Fieldwork for HSE 2011: January 2011 - March 2012
sink(file = paste0(output.dir(), "log.txt"),
     append = T, 
     type = "output",
     split = F)

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


sink()
gc()
