#cmpfile("./output.R")
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

out <- prepare_outcomes(POP)

# collect risk factors
RF <- collect_RF(POP)
nam <- paste0(init.year:(init.year + yearstoproject - 1))
fwrite_safe(summary_RF_fn(RF, nam, "No_HC"), output.dir("sc_RF.csv"))
fwrite_safe(summary_RF_fn(RF, nam, "No_HC", "qimd"), output.dir("sc_RF_qimd.csv"))
fwrite_safe(summary_RF_fn(RF, nam, "No_HC", "sex"), output.dir("sc_RF_sex.csv"))
fwrite_safe(summary_RF_fn(RF, nam, "No_HC", "age"), output.dir("sc_RF_age.csv"))

# shape output
out <- mask_age(out)

out_summary      <- summary_fn(out, nam, "No_HC")
out_summary_qimd <- summary_fn(out, nam, "No_HC", "qimd")
out_summary_sex  <- summary_fn(out, nam, "No_HC", "sex")
out_summary_age  <- summary_fn(out, nam, "No_HC", "age")

fwrite_safe(out_summary,      output.dir("No_HC_results.csv"))
fwrite_safe(out_summary_qimd, output.dir("No_HC_results_qimd.csv"))
fwrite_safe(out_summary_sex,  output.dir("No_HC_results_sex.csv"))
fwrite_safe(out_summary_age,  output.dir("No_HC_results_age.csv"))

rm(RF, out, nam)

