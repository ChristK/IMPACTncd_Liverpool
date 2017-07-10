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
cat("Summarising output...\n")
if (!exists("clusternumber")) clusternumber <- 4L
clusternumber <- ifelse(Sys.info()[1] == "Windows", 1L, clusternumber)
plan(multiprocess)
options(future.globals.maxSize = 32 * 1024^3)
pr <- c(0.5, 0.025, 0.975, .2, 0.8, 0.25, 0.75) # percentiles for summaries

remove_faulty_iterations(output_dir("Intermediate"))

dir.create(
  path = output_dir("Health_Economics/"),
  recursive = T, 
  showWarnings = F
)

dir.create(
  path = output_dir("Graphs/"),
  recursive = T, 
  showWarnings = F
)

dir.create(
  path = output_dir("RF/"),
  recursive = T,
  showWarnings = F
)

dir.create(
  path = output_dir("Validation/"),
  recursive = T,
  showWarnings = F
)

# Health econ ------------------------------------------
all.files <- as.list(
  dir(
    path = output_dir("Intermediate"), 
    pattern = glob2rx("*_sc*_results.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

process_results(results, 3.5, NULL)

all.files <- as.list(
  dir(
    path = output_dir("Intermediate/"), 
    pattern = glob2rx("*_sc*_results_qimd.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results_qimd <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

process_results(results_qimd, 3.5, "qimd")

# rigit scores
results_qimd[, cumul.population := cumsum(population/pop.fraction), keyby = .(scenario, mc, qimd)]
setkey(results_qimd, mc, qimd, scenario)
results_qimd[, qimd2 := cumsum(cumul.population)/sum(cumul.population),
             by = .(mc, year, scenario)]
results_qimd[, qimd2 := (qimd2 + c(0, qimd2[seq_len(.N - 1)]))/2,
             by = .(mc, year, scenario)]
results_qimd[, qimd2.median := median(qimd2, na.rm = T), by = .(year, scenario, qimd)]

# Absolute equity slope index
f_sei <- future({
  results_qimd[, .(sei = fastLm(
  I(cumul.net.utility.disc/pop.fraction) ~ qimd2
)$coefficients["qimd2"]),
by = .(mc, scenario, year)]
})

# Relative equity slope index (exclude qimd = 1)
results_qimd[, proportional.utility.disc := cumul.net.utility.disc/cumul.qalys.disc]
f_rei <- future({
  results_qimd[qimd != "1", .(rei = fastLm(
  proportional.utility.disc ~ qimd2
)$coefficients["qimd2"]),
by = .(mc, scenario, year)]
})

results_qimd[, c("cumul.population", "qimd2.median", "proportional.utility.disc") := NULL]

results_qimd[value(f_sei), on = c("year", "scenario", "mc"), `:=` (sei = i.sei)]
results_qimd[value(f_rei), on = c("year", "scenario", "mc"), `:=` (rei = i.rei)]

results[results_qimd, on = c("year", "scenario", "mc"), `:=` (rei = rei, sei = sei)]

fwrite_safe(results, output_dir("results.csv"), clear_intermediate_outputs)
fwrite_safe(results_qimd, output_dir("results_qimd.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(results, pr, pop.fraction, 0), output_dir("Health_Economics/summary.csv"))
fwrite_safe(summarise_results(results_qimd, pr, pop.fraction, 0, "qimd"), output_dir("Health_Economics/summary_qimd.csv"))

all.files <- as.list(
  dir(
    path = output_dir("Intermediate/"), 
    pattern = glob2rx("*_sc*_results_sex.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results_sex <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)
process_results(results_sex, 3.5, "sex")
fwrite_safe(results_sex, output_dir("results_sex.csv"), clear_intermediate_outputs)

fwrite_safe(summarise_results(results_sex, pr, pop.fraction, 0, "sex"), output_dir("Health_Economics/summary_sex.csv"))

all.files <- as.list(
  dir(
    path = output_dir("Intermediate/"), 
    pattern = glob2rx("*_sc*_results_age.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results_age <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

process_results(results_age, 3.5, "agegroup")

fwrite_safe(results_age, output_dir("results_age.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(results_age, pr, pop.fraction, 0, "agegroup"), output_dir("Health_Economics/summary_age.csv"))

# Graphs -------------------------------------------------------------------
study.year <- (init.year + yearstoproject - 1L)
results[, scenario := factor_(scenario)]
future({
  gg <- 
  ggplot(results[year == study.year],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2, show.legend = F) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median, y = cumul.net.cost.disc_median, 
                 col = scenario),
             data = unique(results[year == study.year,
                                   .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                     cumul.net.cost.disc_median = median(cumul.net.cost.disc, na.rm = T)/pop.fraction), 
                                   by = scenario]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median, y = cumul.net.cost.disc_median, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year,
                                        .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                          cumul.net.cost.disc_median = median(cumul.net.cost.disc, na.rm = T)/pop.fraction), 
                                        by = scenario]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)") +
  scale_y_continuous(name="Incremental cumulative costs (£)", labels = gbp ) +
  scale_colour_discrete(drop=TRUE, breaks = levels(results$scenario), limits = levels(results$scenario)) +
  my_theme 

ggsave("plane_ce.pdf", gg, path = output_dir("Graphs"), dpi = 600, scale = 0.7, width = 8, height = 5)
})

# probability of a policy beeing cost saving
future({
  gg <- 
  ggplot(results[, .(year, prop_if(cumul.net.cost.disc<=0)), 
                 by = .(scenario, year)],
         aes(x = year, y = V2, col = scenario)) +
  geom_point(size = 1, alpha = 5/5) +
  geom_line(size = 1, alpha = 3/5) +
  geom_hline(aes(yintercept = 0.8), show.legend = F) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Probability of cost saving policy", labels = percent ) +
  scale_colour_discrete(drop=TRUE, breaks = levels(results$scenario), limits = levels(results$scenario)) + my_theme 

ggsave("cost_sav_prob.pdf", gg, path = output_dir("Graphs"), dpi = 600, scale = 0.7, width = 8,
       height = 5)
})

# probability of a policy beeing cost effective at $50000 per QALY
future({
  gg <- 
  ggplot(results[, .(year, prop_if(cumul.nmb.disc > 0)), by = .(scenario, year)],
         aes(x = year, y = V2, col = scenario)) +
  geom_point(size = 1, alpha = 5/5) +
  geom_line(size = 1, alpha = 3/5) +
  geom_hline(aes(yintercept = 0.8), show.legend = F) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Probability of cost-effective policy", labels = percent ) +
  scale_colour_discrete(drop=TRUE, breaks = levels(results$scenario), limits = levels(results$scenario)) + my_theme 

ggsave("cost_eff_prob.pdf", gg, path = output_dir("Graphs"), dpi = 600, scale = 0.7, width = 8,
       height = 5) 
})
# Net monetary benefit for 100000$ per qaly
future({
  gg <- 
  ggplot(results[, .(year, quantile(cumul.nmb.disc, c(0.5)) / pop.fraction,
                   quantile(cumul.nmb.disc, c(0.025))/pop.fraction,
                   quantile(cumul.nmb.disc, c(0.975))/pop.fraction),
                 by = .(scenario, year)],
         aes(x = year, y = V2, ymin = V3, ymax = V4, col = scenario, fill = scenario)) +
  geom_point(size = 1, alpha = 5/5, show.legend = F) +
  geom_line(size = 1, alpha = 5/5) +
  geom_ribbon(alpha = 1/5, linetype = 0, show.legend = F) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Net monetary benefit", labels = gbp ) +
  scale_colour_discrete(drop=TRUE, breaks = levels(results$scenario), limits = levels(results$scenario)) + my_theme 

ggsave("nmb.pdf", gg, path = output_dir("Graphs"), dpi = 600, scale = 0.7, width = 8,
       height = 5) 
})
# relative equity
future({
gg <- 
  ggplot(results[year == study.year ],
         aes(x = cumul.net.utility.disc/pop.fraction, y = rei, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median, y = rei.median, 
                 col = scenario),
             data = unique(results[year == study.year,
                                   .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                     rei.median = median(rei, na.rm = T)),
                                   by = scenario]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median, y = rei.median, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year,
                                        .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                          rei.median = median(rei, na.rm = T)),
                                        by = scenario]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)") +
  scale_y_continuous(name="Relative equity index") +
  scale_colour_discrete(drop=TRUE, breaks = levels(results$scenario), limits = levels(results$scenario)) +
  my_theme 
ggsave("rei.pdf", gg, path = output_dir("Graphs"), dpi = 600, scale = 0.7, width = 8,
       height = 5)
})

# absolute equity
future({
gg <- 
  ggplot(results[year == study.year ],
         aes(x = cumul.net.utility.disc/pop.fraction, y = sei, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median, y = sei.median, 
                 col = scenario),
             data = unique(results[year == study.year,
                                   .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                     sei.median = median(sei, na.rm = T)),
                                   by = scenario]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median, y = sei.median, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year,
                                        .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                          sei.median = median(sei, na.rm = T)),
                                        by = scenario]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)") +
  scale_y_continuous(name="Absolute equity index") +
  scale_colour_discrete(drop=TRUE, breaks = levels(results$scenario), limits = levels(results$scenario)) +
  my_theme 

ggsave("sei.pdf", gg, path = output_dir("Graphs"), dpi = 600, scale = 0.7, width = 8,
       height = 5)
})

# RF ----------------------------------------------------------------------
all.files <- as.list(
  dir(
    path = output_dir("Intermediate"), 
    pattern = glob2rx("*_RF.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

fwrite_safe(results, output_dir("RF.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(results, pr, 1, 3, NULL), output_dir("RF/RF_summary.csv"))

all.files <- as.list(
  dir(
    path = output_dir("Intermediate"), 
    pattern = glob2rx("*_RF_sex.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results_sex <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

fwrite_safe(results_sex, output_dir("RF_sex.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(results_sex, pr, 1, 3, "sex"), output_dir("RF/RF_summary_sex.csv"))

all.files <- as.list(
  dir(
    path = output_dir("Intermediate"), 
    pattern = glob2rx("*_RF_qimd.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results_qimd <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

fwrite_safe(results_qimd, output_dir("RF_qimd.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(results_qimd, pr, 1, 3, "qimd"), output_dir("RF/RF_summary_qimd.csv"))

all.files <- as.list(
  dir(
    path = output_dir("Intermediate"), 
    pattern = glob2rx("*_RF_age.csv"), 
    full.names = T, 
    recursive = F
  )
) 

results_age <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)


fwrite_safe(results_age, output_dir("RF_age.csv"), clear_intermediate_outputs)
# fwrite_safe(summarise_RF(results_age, pr, 1, 3,, "age"), output_dir("RF/RF_summary_age.csv"))

# Validation --------------------------------------------------------------
all.files <- as.list(
  dir(
    path = output_dir("Intermediate/"), 
    pattern = glob2rx("*_No_HC_results.csv"), 
    full.names = T, 
    recursive = F
  )
) 

validat <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

fwrite_safe(validat, output_dir("Validation/results.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(validat, pr, pop.fraction, 0, NULL), output_dir("Validation/summary.csv"))

validat_sum <- 
  validat[between(year, 2011, 2015), c(output_quants(chd.deaths, "chd.deaths"),
              output_quants(stroke.deaths, "stroke.deaths"),
              output_quants(other.deaths + chd.deaths + stroke.deaths,
                            "all.deaths")),
          keyby = .(year)]
validat_sum <- 
  melt(validat_sum, 1, variable.name = "type", value.name = "deaths")
validat_sum[, source := "IMPACTncd"]

observed_mort <- fread("./Validation/mortality_observed.csv", header = T)
observed_mort <- 
  melt(observed_mort, 1:4, variable.name = "qimd", value.name = "deaths")
tt <- observed_mort[, .(deaths = sum(deaths)), by = .(year, type)]
tt[, source := "Observed"]
validat_sum <- rbind(validat_sum, tt)

gg <- qplot(year, deaths, colour = source,
            data = validat_sum[type == "chd.deaths_50.0%"],
      main = "CHD mortality")
ggsave(output_dir("Validation/chd_mort_valid.pdf"), gg)
gg <- qplot(year, deaths, colour = source,
            data = validat_sum[type == "stroke.deaths_50.0%"],
      main = "Stroke mortality")
ggsave(output_dir("Validation/stroke_mort_valid.pdf"), gg)
gg <- qplot(year, deaths, colour = source,
            data = validat_sum[type == "all.deaths_50.0%"],
      main = "All-cause mortality")
ggsave(output_dir("Validation/all_mort_valid.pdf"), gg)

all.files <- as.list(
  dir(
    path = output_dir("Intermediate/"), 
    pattern = glob2rx("*_No_HC_results_qimd.csv"), 
    full.names = T, 
    recursive = F
  )
) 

validat_qimd <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

fwrite_safe(validat_qimd, output_dir("/Validation/results_qimd.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(validat_qimd, pr, pop.fraction, 0, "qimd"), output_dir("Validation/summary_qimd.csv"))

validat_sum_qimd <- 
  validat_qimd[between(year, 2011, 2015), c(output_quants(chd.deaths, "chd.deaths"),
                                       output_quants(stroke.deaths, "stroke.deaths"),
                                       output_quants(other.deaths + chd.deaths + stroke.deaths,
                                                     "all.deaths")),
          keyby = .(year, qimd)]
validat_sum_qimd <- 
  melt(validat_sum_qimd, 1:2, variable.name = "type", value.name = "deaths")
validat_sum_qimd[, source := "IMPACTncd"]

observed_mort <- fread("./Validation/mortality_observed.csv", header = T)
observed_mort <- 
  melt(observed_mort, 1:4, variable.name = "qimd", value.name = "deaths")
tt <- observed_mort[, .(deaths = sum(deaths)), by = .(year, qimd, type)]
tt[, source := "Observed"]
validat_sum_qimd <- rbind(validat_sum_qimd, tt)

gg <- qplot(year, deaths, colour = source, facets = .~qimd,
            data = validat_sum_qimd[type == "chd.deaths_50.0%"],
            main = "CHD mortality")
ggsave(output_dir("Validation/chd_mort_qimd_valid.pdf"), gg)
gg <- qplot(year, deaths, colour = source, facets = .~qimd,
            data = validat_sum_qimd[type == "stroke.deaths_50.0%"],
            main = "Stroke mortality")
ggsave(output_dir("Validation/stroke_mort_qimd_valid.pdf"), gg)
gg <- qplot(year, deaths, colour = source, facets = .~qimd,
            data = validat_sum_qimd[type == "all.deaths_50.0%"],
            main = "All-cause mortality")
ggsave(output_dir("Validation/all_mort_qimd_valid.pdf"), gg)

all.files <- as.list(
  dir(
    path = output_dir("Intermediate/"), 
    pattern = glob2rx("*_No_HC_results_sex.csv"), 
    full.names = T, 
    recursive = F
  )
) 

validat_sex <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)

fwrite_safe(validat_sex, output_dir("/Validation/results_sex.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(validat_sex, pr, pop.fraction, 0, "sex"), output_dir("Validation/summary_sex.csv"))

all.files <- as.list(
  dir(
    path = output_dir("Intermediate/"), 
    pattern = glob2rx("*_No_HC_results_age.csv"), 
    full.names = T, 
    recursive = F
  )
) 

validat_age <- rbindlist(
  mclapply(
    all.files,
    fread,
    mc.cores = clusternumber
  ),
  T, T
)


fwrite_safe(validat_age, output_dir("Validation/results_age.csv"), clear_intermediate_outputs)
fwrite_safe(summarise_results(validat_age, pr, pop.fraction, 0, "agegroup"), output_dir("Validation/summary_age.csv"))

# Clear intermediate files ------------------------------------------------
if (clear_intermediate_outputs == T) {
  unlink(output_dir("Intermediate/"), recursive = T, force = T)
}
# rm(
#   list = setdiff(
#     ls(),
#     lsf.str()
#   )
# ) # remove everything but functions


