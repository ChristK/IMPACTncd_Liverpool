#cmpfile("./2dmc.R")
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


cat("Sample RR values for 2d Monte Carlo\n\n")
# coefficients for salt model from the MC simulation
cvd.lag.mc <- cvd.lag.l[[iterations]]
#cancer.lag <- cancer.lag.l[[iterations]]

load(file="./Lagtimes/salt.rq.coef.rda")
salt.rq$coefficients <- salt.rq.coef[[iterations]]
#salt.rq$coefficients <- apply(simplify2array(salt.rq.coef), 1:2, mean) # mean of MC

atorv.eff.mc          <- atorv.eff.l[[iterations]]
persistence.mc        <- persistence.l[[iterations]]
adherence.mc          <- adherence.l[[iterations]]
ssb.tax.effect.mc     <- ssb.tax.effect.l[.id == iterations]
salt.reform.effect.mc <- salt.reform.effect.l[.id == iterations]
smok.tcs.effect.mc    <- smok.tcs.effect.l[.id == iterations]

if ("CHD" %in% diseasestoexclude) {
  chd.cost.death.mc <- chd.cost.death.l[[iterations]]
  tobacco.rr.chd.mc <- chd.tobacco.rr.l[.id == iterations]
  chd.ets.rr.mc     <- chd.ets.rr.l[[iterations]]
  sbp.rr.chd.mc     <- chd.sbp.rr.l[.id == iterations]
  chol.rr.chd.mc    <- chd.chol.rr.l[.id == iterations]
  chd.bmi.rr.mc     <- chd.bmi.rr.l[.id == iterations]
  chd.diab.rr.mc    <- chd.diab.rr.l[.id == iterations]
  chd.fv.rr.mc      <- chd.fv.rr.l[[iterations]]
  pa.rr.chd.mc      <- chd.pa.rr.l[.id == iterations]
}

if ("stroke" %in% diseasestoexclude) {
  stroke.cost.death.mc <- stroke.cost.death.l[[iterations]]
  tobacco.rr.stroke.mc <- stroke.tobacco.rr.l[.id == iterations]
  stroke.ets.rr.mc     <- stroke.ets.rr.l[[iterations]]
  sbp.rr.stroke.mc     <- stroke.sbp.rr.l[.id == iterations]
  chol.rr.stroke.mc    <- stroke.chol.rr.l[.id == iterations]
  stroke.bmi.rr.mc     <- stroke.bmi.rr.l[.id == iterations]
  stroke.diab.rr.mc    <- stroke.diab.rr.l[.id == iterations]
  stroke.fv.rr.mc      <- stroke.fv.rr.l[[iterations]]
  pa.rr.stroke.mc      <- stroke.pa.rr.l[.id == iterations]
}


# Health economics --------------------------------------------------------
utility.mc <- utility.l[.id == iterations]
cost.mc <- cost.l[.id == iterations]
