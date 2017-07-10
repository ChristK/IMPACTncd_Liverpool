
# Preample ---------------------------------------------------
clusternumber <- 10L
 
get_dropbox_path <- function(x = character(0)) {
  # if (!require(jsonlite)) stop("You need to install jsonlite package.")
  if (!require(jsonlite)) install.packages("jsonlite")
  if (Sys.info()[1] == "Windows")
  {
    if (file.exists(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json")))
    {
      dropbox_path <- read_json(paste0(Sys.getenv("APPDATA"), "/Dropbox/info.json"))$personal$path
    } else dropbox_path <- read_json(paste0(Sys.getenv("LOCALAPPDATA"), "/Dropbox/info.json"))$personal$path
  } else {
    dropbox_path <- read_json("~/.dropbox/info.json")$personal$path
  }
  return(paste0(dropbox_path, "/", x))
}

if (grepl("incline", Sys.info()["nodename"], fixed = T))
{
  setwd("~/IMPACTncd Liverpool/")
  clusternumber <- ifelse(clusternumber > 30, 30, clusternumber)
  output_dir <- function(x = character(0))
    (paste0("./Output/", x))
  synthpop_dir <- "./SynthPop"
} else if (Sys.info()["nodename"] == "dh012142") {
  setwd(get_dropbox_path("Models/IMPACTncd Liverpool/"))
  clusternumber <- ifelse(clusternumber > 25, 25, clusternumber)
  output_dir <-
    function(x = character(0))
      (paste0("/mnt/storage_slow/Model_Outputs/Liverpool_Model/", x))
  synthpop_dir <-
    "~/Synthetic Populations/Northwest England 2011/"
} else {
  setwd(get_dropbox_path("Models/IMPACTncd Liverpool/"))
  output_dir <- function(x = character(0))
    (paste0("./Output/", x))
  synthpop_dir <- "./SynthPop"
}


require(Cairo)
if (Sys.info()[1] == "Windows") {
  Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.20/bin/gswin64c.exe")
  require(extrafont)
  # font_import() # to be run only once in each system
  #loadfonts(device = "win", quiet = T)
  loadfonts()
  #fonttable()
} 
require(data.table)
require(ggplot2)
require(ggrepel)
require(scales)
require(ggthemes)
require(Rcpp)
require(RColorBrewer)
require(tikzDevice)
options(tikzDocumentDeclaration = "\\documentclass[tikz]{standalone}")
options(tikzMetricPackages = c(
  "\\usepackage[utf8]{inputenc}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage{siunitx}",
  "\\usetikzlibrary{calc}",
  "\\usepackage{amssymb}",
  "\\sisetup{mode = text, detect-all, per-mode=symbol}"))
options("tikzLatexPackages" =
          c(getOption("tikzLatexPackages"),"\\usepackage{siunitx}"))

# output quantiles
output_quants <- 
    function(x, nam = get_object_name(x), pop_extrapol = pop.fraction, ...) {
      pr <- c(0.025, .2, 0.25, 0.5, 0.75, 0.8, 0.975)
      pr_nam <- percent(pr) #Very Slow
      out <- transpose(list(round(quantile(x, pr, T,
                                           F) / pop_extrapol, ...)))
      setattr(out, "names", paste0(nam, "_", pr_nam))
      return(out)
    }

output_quants_lbl <-
  function(x,
           nam = get_object_name(x),
           pop_extrapol = pop.fraction,
           prob = c(0.5, 0.025, 0.975),
           digits = 2,
           ...) {
    out <- as.list(signif(quantile(x, prob,
                                   na.rm = T) / pop_extrapol, digits = digits))
    out <-
      list(paste0(
        prettyNum(out[1], ","),
        " (",
        prettyNum(out[2], ","),
        " to ",
        prettyNum(out[3], ","),
        ")"
      ))
    #names(out) <- paste0(nam)
    setattr(out, "names", paste0(nam))
    return(out)
  }


qimd_labeller <- labeller(qimd = c("1" = "QIMD 1\nleast deprived",
                                   "2" = "QIMD 2",
                                   "3" = "QIMD 3",
                                   "4" = "QIMD 4",
                                   "5" = "QIMD 5\nmost deprived"))
gbp <- dollar_format(prefix = "£", suffix = "")

# normalise a vector to 0,1 range
sourceCpp("functions.cpp", cacheDir = "./.CppCache/")
identical_elements <- 
    function(x, tol = .Machine$double.eps ^ 0.5) {
      stopifnot(is.numeric(x))
      fequal(x, tol)
    }

normalise <-
    function(x, ...){
      stopifnot(is.numeric(x))
      if (identical_elements(x)) return(1)
      else return(fnormalise(x))
    }

my_theme <-   theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = rel(1), vjust = 0),
        axis.title = element_text(size = rel(1)),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(size = rel(0.8), vjust = 0.5, angle = 90),
        axis.text.y = element_text(size = rel(0.8), vjust = 0),
        strip.text.x = element_text(rel(.8)),
        strip.text.y = element_text(rel(.8)),
        strip.background = element_rect(colour = "lavenderblush4", fill = "lavenderblush2"),
        legend.title = element_blank(), # no legent title
        legend.position = "right",
        #legend.justification = c(0,1), legend.position = c(0,1),
        legend.key = element_rect(size = 1, colour = "transparent", fill = "transparent"),
        legend.key.size = unit(1, "lines")) 

# Export tikz
save_tikz <-
  function(filename,
           graph,
           path,
           width = 11.690 / 1.4,
           height = 8.270 / 1.4) {
    tikz(
      file = paste0(path, filename),
      width = width,
      height = height,
      lwdUnit = 72.27 / 96,
      standAlone = T,
      sanitize = T,
      packages = c(
        "\\usepackage{tikz}",
        "\\usepackage[T1]{fontenc}",
        "\\usepackage[active,tightpage,psfixbb]{preview}",
        "\\PreviewEnvironment{pgfpicture}",
        "\\setlength\\PreviewBorder{0pt}"
      )
    )
    print(graph)
    dev.off()
  }

# Load data
results <- fread(output_dir("results.csv"))
if (!exists("pop.fraction"))
{
  tt <- readLines(output_dir("simulation parameters.txt"))
  pop.fraction <-
    as.numeric(substring(tt[[grep(glob2rx("Population fraction = *"), tt)]], 23))
}

# Process data ----------------------------------------------------
results[scenario == "sc0", scenario := "Current"]
results[scenario == "sc1", scenario := "Uptake"]
results[scenario == "sc2", scenario := "Coverage"]
results[scenario == "sc5", scenario := "Concentrated"]
results[scenario == "sc6a", scenario := "Prescription"]
results[scenario == "sc7a", scenario := "Lifestyle"]
results[scenario == "sc7b", scenario := "Optimal"]
results[scenario == "sc11", scenario := "Maximising TCS"]
results[, scenario := factor(scenario)]

results[, cumul.net.utility.disc.median := median(cumul.net.utility.disc, na.rm = T), by = .(year, scenario)]
results[, cumul.net.cost.disc.median := median(cumul.net.cost.disc, na.rm = T), by = .(year, scenario)]

# Graphs ------------------------------------------------------------------

study.year <- 2031 # 15 years post intervention
study.sc <- "Current"
colourCount <- length(unique(levels(results$scenario)))
getPalette <- colorRampPalette(brewer.pal(9, "Set1")) # create longer colour palette by interpolation
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2, show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 3, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% study.sc,
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme + theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_current15_noSE.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2, show.legend = F) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 3, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                       col = scenario, label = scenario),
                   data = unique(results[year == study.year & scenario %in% study.sc,
                                         .(cumul.net.utility.disc.median, 
                                           cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme  + theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0)) 

ggsave("CEplane_current15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

study.sc <- c("Current", "Uptake")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme + theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_uptake15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")


study.sc <- c("Current", "Uptake", "Coverage")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_coverage15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

study.sc <- c("Current", "Uptake", "Coverage", "Concentrated")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = 450, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_concetrated15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")


study.sc <- c("Current", "Uptake", "Coverage", "Concentrated", "Prescription")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_prescription15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

study.sc <- c("Current", "Uptake", "Coverage", "Concentrated", "Prescription", "Lifestyle")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_lifestyle15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

study.sc <- c("Current", "Optimal")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% study.sc,
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = 850, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 2000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-1.5e7, 1.3e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_optimal15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

results[scenario == "sc8", scenario := "???"]
results[, scenario := factor(scenario)]

study.sc <- c("Current", "Optimal", "???")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = 50, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 2000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-1.5e7, 1.3e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_salt15_q.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

results[scenario == "???", scenario := "Salt"]
results[, scenario := factor(scenario)]

study.sc <- c("Current", "Optimal", "Salt")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = 50, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 2000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-1.5e7, 1.3e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_salt15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

results[scenario == "sc12a", scenario := "???"]
results[, scenario := factor(scenario)]

study.sc <- c("Current", "Optimal", "Salt", "???")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = 50, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 6000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-3.5e7, 1.3e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_diet15_q.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

results[scenario == "???", scenario := "Diet"]
results[, scenario := factor(scenario)]

study.sc <- c("Current", "Optimal", "Salt", "Diet")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = 50, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 6000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-3.5e7, 1.3e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_diet15.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

gg <- 
  ggplot(results[year == study.year & scenario %in% c("Current")],
         aes(x = cumul.net.utility.disc, y = cumul.net.cost.disc, col = scenario)) +
  geom_point(alpha = 0) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 1000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-2e6, 1e7), labels = gbp ) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))
        #legend.position = "") 


ggsave("CEplane.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

study.sc <- c("Current", "Maximising TCS")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = cumul.net.cost.disc/pop.fraction, col = scenario)) +
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) +
  geom_abline(slope = 2e4, intercept = 0, linetype = 2) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median, 
                                     cumul.net.cost.disc.median, scenario)]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median/pop.fraction, y = cumul.net.cost.disc.median/pop.fraction, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% last(study.sc),
                                        .(cumul.net.utility.disc.median, 
                                          cumul.net.cost.disc.median, scenario)]), inherit.aes = F,
                  nudge_x = 50, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)", limits = c(-500, 6000)) +
  scale_y_continuous(name="Incremental cumulative costs (£)", limits = c(-3.5e7, 1.3e7), labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("CEplane_TCS15_noSE.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

# Graphs -------------------------------------------------------------
# absolute equity
study.sc <- c("Current", "Optimal", "Diet")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = sei, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median, y = sei.median, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                     sei.median = median(sei, na.rm = T)),
                                   by = scenario]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median, y = sei.median, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% study.sc,
                                        .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                          sei.median = median(sei, na.rm = T)),
                                        by = scenario]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)") +
  scale_y_continuous(name="Absolute equity index") +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("sei.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

# relative equity
study.sc <- c("Current", "Optimal", "Diet")
gg <- 
  ggplot(results[year == study.year & scenario %in% study.sc],
         aes(x = cumul.net.utility.disc/pop.fraction, y = rei, col = scenario)) +
  geom_hline(aes(yintercept = 0), show.legend = F) + geom_vline(aes(xintercept = 0), show.legend = F) +
  stat_ellipse(type = "norm", show.legend = F) +
  geom_point(aes(x = cumul.net.utility.disc.median, y = rei.median, 
                 col = scenario),
             data = unique(results[year == study.year & scenario %in% study.sc,
                                   .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                     rei.median = median(rei, na.rm = T)),
                                   by = scenario]),
             size = 2, alpha = 5/5,
             inherit.aes = F) +
  geom_text_repel(aes(x = cumul.net.utility.disc.median, y = rei.median, 
                      col = scenario, label = scenario),
                  data = unique(results[year == study.year & scenario %in% study.sc,
                                        .(cumul.net.utility.disc.median = median(cumul.net.utility.disc, na.rm = T)/pop.fraction, 
                                          rei.median = median(rei, na.rm = T)),
                                        by = scenario]), inherit.aes = F,
                  nudge_x = -400, show.legend = F) +
  scale_x_continuous(name="Incremental cumulative effects (QALYs)") +
  scale_y_continuous(name="Relative equity index", labels = percent) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("rei.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

# probability of a policy being cost saving
study.sc <- c("Current", "Optimal", "Diet")
gg <- 
  ggplot(results[scenario %in% study.sc, .(year, prop_if(cumul.net.cost.disc<=0)), 
                 by = .(scenario, year)],
         aes(x = year, y = V2, col = scenario)) +
  geom_point(size = 1, alpha = 5/5) +
  geom_line(size = 1, alpha = 3/5) +
  geom_hline(aes(yintercept = 0.8), show.legend = F) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Probability of cost saving policy", labels = percent ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("cost_sav_prob.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")


# probability of a policy being cost effective at £20000 per QALY
study.sc <- c("Current", "Optimal", "Diet")
gg <- 
  ggplot(results[scenario %in% study.sc, .(year, prop_if(cumul.nmb.disc > 0)), by = .(scenario, year)],
         aes(x = year, y = V2, col = scenario)) +
  geom_point(size = 1, alpha = 5/5) +
  geom_line(size = 1, alpha = 3/5) +
  geom_hline(aes(yintercept = 0.8), show.legend = F) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Probability of cost-effective policy", labels = percent ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("cost_eff_prob.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")

# Net monetary benefit for 100000$ per qaly
gg <- 
  ggplot(results[scenario %in% study.sc, .(year, quantile(cumul.nmb.disc, c(0.5)) / pop.fraction,
                     quantile(cumul.nmb.disc, c(0.025))/pop.fraction,
                     quantile(cumul.nmb.disc, c(0.975))/pop.fraction),
                 by = .(scenario, year)],
         aes(x = year, y = V2, ymin = V3, ymax = V4, col = scenario, fill = scenario)) +
  geom_point(size = 1, alpha = 5/5, show.legend = F) +
  geom_line(size = 1, alpha = 5/5) +
  geom_ribbon(alpha = 1/5, linetype = 0, show.legend = F) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Net monetary benefit", labels = gbp ) +
  scale_colour_manual(values = getPalette(colourCount), drop=TRUE, breaks = study.sc, limits = levels(results$scenario)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0))

ggsave("nmb.tiff", gg, path = sav.dir, dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")


# relative equity
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
ggsave("rei.pdf", gg, path = "./Output/Graphs/", dpi = 600, scale = 0.7, width = 8,
       height = 5)


# Experimental graphs -----------------------------------------------------
#library("shinyjs")
study.year <- 2031
study.sc <- c("Current", "Uptake", "Coverage", "Concentrated", "Prescription", "Lifestyle", "Optimal", "Diet")

tt <- results[year == study.year & scenario %in% study.sc,
              .("Cost\nsaving" = sum(cumul.net.cost.disc < 0)/.N, # inverse <0 because higher cost is bad
                "QALYs    " = sum(cumul.net.utility.disc > 0)/.N,
                "Cost\neffective" = sum(cumul.nmb.disc > 0)/.N # probability of cheaper than £20000 per qualy
              ),
              keyby= .(scenario)]

ttt <- results[year == study.year & scenario %in% study.sc, #any qimd other than 1 will do
              .("     Absolute\n     equity" = sum(sei > 0)/.N, # 
                "Relative\nequity" = sum(rei > 0)/.N
                ),
              keyby= .(scenario)]

tt <- melt(tt, 1)
ttt <- melt(ttt, 1)
tt <- rbind(tt, ttt)
tt[, scenario := factor(scenario, levels = study.sc)]
tt[, variable := factor(variable)]

gg <- 
  ggplot(tt, aes(x = variable, y = value)) + 
  geom_col(fill = "#1E90FF", alpha = 0.5,  
           width = 1, show.legend = FALSE, color = "white") +
  geom_hline(yintercept = 0.8, 
             colour = "red", size = 1, lty = 1) +
  geom_hline(yintercept = c(0.2, 0.4, 0.6, 1), 
             colour = "#949494", size = 1, lty = 3) +
  geom_vline(xintercept = seq(0.5, 5.5, 1), # as many as variables
             colour = "#949494", size = 1, lty = 1) +
  facet_wrap(~scenario, ncol = 4) +
  coord_polar() + 
  #scale_fill_manual(values = c("#e4e4e4", "#00a9e0")) + 
  scale_y_continuous(
    limits = c(0, 1), 
    breaks = c(0.2, 0.4, 0.6, 0.8, 1)) + 
  labs(x = "", y = "") +

  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = rel(1), vjust = 0),
        axis.title = element_text(size = rel(1)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2)),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 18),
        strip.background = element_rect(colour = "lavenderblush4", fill = "lavenderblush2"),
        legend.title = element_blank(), # no legent title
        legend.position = "none",
        panel.spacing = unit(0.5, "lines")) 


ggsave("spider.tiff", gg, path = "./Output/Graphs/", dpi = 600, scale = 2, width = 8, height = 6,
       compression = "lzw", antialias = "subpixel", type="cairo")

  
tt <- results[ scenario %in% study.sc,
              .("Cost saving" = sum(cumul.net.cost.disc < 0)/.N, # inverse <0 because higher cost is bad
                "QALYs" = sum(cumul.net.utility.disc > 0)/.N,
                "Cost effective" = sum(cumul.nmb.disc > 0)/.N # probability of cheaper than £20000 per qualy
              ),
              keyby= .(scenario, year)]

ttt <- results[ scenario %in% study.sc , #any qimd other than 1 will do
                    .("Absolute equity" = sum(sei > 0)/.N, # 
                      "Relative equity" = sum(rei > 0)/.N
                    ),
                    keyby= .(scenario, year)]

tt <- melt(tt, 1:2)
ttt <- melt(ttt, 1:2)
tt <- rbind(tt, ttt)
tt[, scenario := factor(scenario, levels = study.sc)]
tt[, variable := factor(variable)]
tt[, value_norm := value - 0.8, by = .(variable)]
tt[value_norm < 0 , value_norm := -normalise(-value_norm)]
tt[value_norm > 0 , value_norm := normalise(value_norm)]
setkey(tt, year)

gg <- 
  ggplot(tt, aes(x = scenario, y = year, fill = value_norm, colour = value_norm)) + 
  geom_col(alpha = 1,  #position = "fill",
           width = 1, show.legend = FALSE) +
  # geom_hline(yintercept = 0.8, 
  #            colour = "red", size = 1, lty = 1) +
  # geom_hline(yintercept = c(2020, 2030),
  #            colour = "#949494", size = 1, lty = 3) +
  geom_vline(xintercept = seq(0.5, 8.5, 1), # as many as variables
             colour = "#949494", size = 1, lty = 1) +
  scale_colour_gradient2(low = "red", high = "blue") +
  scale_fill_gradient2(low = "red", high = "blue") +
  facet_wrap(~variable, ncol = 3) +
  coord_flip() +
  labs(x = "", y = "") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = rel(1), vjust = 0),
        axis.title = element_text(size = rel(1)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2)),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = rel(2), angle = 0),
        axis.text.x = element_blank(),
        strip.background = element_rect(colour = "lavenderblush4", fill = "lavenderblush2"),
        legend.title = element_blank(), # no legent title
        legend.position = "none",
        panel.spacing = unit(0.5, "lines")) 

ggsave("flame.tiff", gg, path = sav.dir, dpi = 600, scale = 1, width = 8, height = 6,
       compression = "lzw", antialias = "subpixel", type="cairo")

# Validation --------------------------------------------------------------

validat <- fread(output_dir("Validation/results_qimd.csv"))

validat <- 
  validat[between(year, 2011, 2015), c(output_quants(chd.deaths, "chd.deaths"),
                                       output_quants(stroke.deaths, "stroke.deaths"),
                                       output_quants(other.deaths + chd.deaths + stroke.deaths,
                                                     "all.deaths")),
          keyby = .(year, qimd)]
validat <- 
  melt(validat, 1:2, variable.name = "type", value.name = "deaths")
validat[, source := "IMPACTncd"]

observed_mort <- fread("./Validation/mortality_observed.csv", header = T)
observed_mort <- 
  melt(observed_mort, 1:4, variable.name = "qimd", value.name = "deaths")
tt <- observed_mort[, .(deaths = sum(deaths)), by = .(year, qimd, type)]
tt[, source := "Observed"]
validat <- rbind(validat, tt)


gg <- ggplot(validat[type == "chd.deaths_50.0%"], aes(year, deaths, colour = source)) +
  geom_point(size = 3, alpha=3/5) +
  facet_grid(. ~ qimd, labeller = qimd_labeller) +
  ylab("CHD deaths") +
  scale_x_continuous(name="Year", limits = c(2011, 2015), breaks = seq(2011, 2015, 2)) +
  scale_colour_discrete(breaks = c("Observed", "IMPACTncd"),
                        labels = c("Observed   ", expression("IMPACT"["NCD"]))) + 
  my_theme + theme(legend.position = "top")

ggsave("chd_mort_qimd_valid.tiff", gg, path = output_dir("Report/Validation/"), dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")
# ggsave("chd_mort_qimd_valid.png", gg, path = "./Output/Validation/", dpi = 600, scale = 0.7, width = 8, height = 5,
#       antialias = "subpixel", type="cairo")
save_tikz("chd_mort_qimd_valid.tex", gg, path = output_dir("Report/Validation/"))

gg <- ggplot(validat[type == "stroke.deaths_50.0%"], aes(year, deaths, colour = source)) +
  geom_point(size = 3, alpha=4/5) +
  facet_grid(. ~ qimd, labeller = qimd_labeller) +
  ylab("Stroke deaths") +
  scale_x_continuous(name="Year", limits = c(2011, 2015), breaks = seq(2011, 2015, 2)) +
  scale_colour_discrete(breaks = c("Observed", "IMPACTncd"),
                        labels = c("Observed   ", expression("IMPACT"["NCD"]))) + 
  my_theme + theme(legend.position = "top")

ggsave("stroke_mort_qimd_valid.tiff", gg, path = output_dir("Report/Validation/"), dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")
# ggsave("stroke_mort_qimd_valid.png", gg, path = "./Output/Validation/", dpi = 600, scale = 0.7, width = 8, height = 5,
#        antialias = "subpixel", type="cairo")
save_tikz("stroke_mort_qimd_valid.tex", gg, path = output_dir("Report/Validation/"))

gg <- ggplot(validat[type == "all.deaths_50.0%"], aes(year, deaths, colour = source)) +
  geom_point(size = 3, alpha=4/5) +
  facet_grid(. ~ qimd, labeller = qimd_labeller) +
  ylab("All-cause deaths") +
  scale_x_continuous(name="Year", limits = c(2011, 2015), breaks = seq(2011, 2015, 2)) +
  scale_colour_discrete(breaks = c("Observed", "IMPACTncd"),
                        labels = c("Observed   ", expression("IMPACT"["NCD"]))) + 
  my_theme + theme(legend.position = "top")

ggsave("all_mort_qimd_valid.tiff", gg, path = output_dir("Report/Validation/"), dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")
# ggsave("all_mort_qimd_valid.png", gg, path = "./Output/Validation/", dpi = 600, scale = 0.7, width = 8, height = 5,
#       antialias = "subpixel", type="cairo")
save_tikz("all_mort_qimd_valid.tex", gg, path = output_dir("Report/Validation/"))


# Tables for report -------------------------------------------------------
validat <- fread(output_dir("/Validation/results.csv"))
discount <- (100 - 3.5)/100

validat[, `:=` (
  qalys.disc = qalys * discount ^ (year - 2016),
  disease.cost.disc = disease.cost * discount ^ (year - 2016)
),
by = c("scenario", "mc")]

setkeyv(validat, c("year"))
validat[, `:=` (
  cumul.chd.incid = cumsum(chd.incid),
  cumul.stroke.incid = cumsum(stroke.incid),
  cumul.cvd.incid = cumsum(cvd.incid),
  cumul.chd.deaths = cumsum(chd.deaths),
  cumul.stroke.deaths = cumsum(stroke.deaths),
  cumul.qalys.disc = cumsum(qalys.disc),
  cumul.disease.cost.disc = cumsum(disease.cost.disc)
),
by = c("scenario", "mc")]

validat[, summary(cvd.incid - (cvd.incid - chd.incid) - (cvd.incid - stroke.incid))]
tt <- validat[, c(
  output_quants_lbl(population, "population", prob = c(0.5, 0.025, 0.975), digits = 3),
  output_quants_lbl(cumul.cvd.incid - cumul.stroke.incid, "cumul.chd.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - stroke.incid)/population, "chd.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - stroke.preval)/population, "chd.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.cvd.incid - cumul.chd.incid, "cumul.stroke.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - chd.incid)/population, "stroke.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - chd.preval)/population, "stroke.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.chd.deaths, "cumul.chd.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*chd.deaths/population, "chd.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.stroke.deaths, "cumul.stroke.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*stroke.deaths/population, "stroke.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.qalys.disc/1e6, "cumul.qalys.disc.millions", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.disease.cost.disc/1e6, "cumul.disease.cost.disc.millions", prob = c(0.5, 0.025, 0.975))
),
keyby = c("year")]
tt <- tt[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year")), variable~year), output_dir("Report/baseline_summary_results.csv")) 



validat <- fread(output_dir("Validation/results_sex.csv"))
discount <- (100 - 3.5)/100

validat[, `:=` (
  qalys.disc = qalys * discount ^ (year - 2016),
  disease.cost.disc = disease.cost * discount ^ (year - 2016)
),
by = c("scenario", "mc")]

setkeyv(validat, c("year", "sex"))
validat[, `:=` (
  cumul.chd.incid = cumsum(chd.incid),
  cumul.stroke.incid = cumsum(stroke.incid),
  cumul.cvd.incid = cumsum(cvd.incid),
  cumul.chd.deaths = cumsum(chd.deaths),
  cumul.stroke.deaths = cumsum(stroke.deaths),
  cumul.qalys.disc = cumsum(qalys.disc),
  cumul.disease.cost.disc = cumsum(disease.cost.disc)
),
by = c("scenario", "mc", "sex")]

validat[, summary(cvd.incid - chd.incid - stroke.incid)]
tt_sex <- validat[, c(
  output_quants_lbl(population, "population", prob = c(0.5, 0.025, 0.975), digits = 3),
  output_quants_lbl(cumul.cvd.incid - cumul.stroke.incid, "cumul.chd.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - stroke.incid)/population, "chd.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - stroke.preval)/population, "chd.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.cvd.incid - cumul.chd.incid, "cumul.stroke.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - chd.incid)/population, "stroke.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - chd.preval)/population, "stroke.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.chd.deaths, "cumul.chd.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*chd.deaths/population, "chd.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.stroke.deaths, "cumul.stroke.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*stroke.deaths/population, "stroke.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.qalys.disc/1e6, "cumul.qalys.disc.millions", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.disease.cost.disc/1e6, "cumul.disease.cost.disc.millions", prob = c(0.5, 0.025, 0.975))
),
keyby = c("year", "sex")]
tt <- tt_sex[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year", "sex")), sex+variable~year), output_dir("Report/baseline_summary_results_sex.csv")) 

# table in the appendix
validat <- fread(output_dir("Validation/results_qimd.csv"))
discount <- (100 - 3.5)/100

validat[, `:=` (
  qalys.disc = qalys * discount ^ (year - 2016),
  disease.cost.disc = disease.cost * discount ^ (year - 2016)
),
by = c("scenario", "mc", "qimd")]

setkeyv(validat, c("year", "qimd"))
validat[, `:=` (
  cumul.chd.incid = cumsum(chd.incid),
  cumul.stroke.incid = cumsum(stroke.incid),
  cumul.cvd.incid = cumsum(cvd.incid),
  cumul.chd.deaths = cumsum(chd.deaths),
  cumul.stroke.deaths = cumsum(stroke.deaths),
  cumul.qalys.disc = cumsum(qalys.disc),
  cumul.disease.cost.disc = cumsum(disease.cost.disc)
  ),
  by = c("scenario", "mc", "qimd")]

validat[, summary(cvd.incid - chd.incid - stroke.incid)]
tt_qimd <- validat[, c(
  output_quants_lbl(population, "population", prob = c(0.5, 0.025, 0.975), digits = 3),
  output_quants_lbl(cumul.cvd.incid - cumul.stroke.incid, "cumul.chd.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - stroke.incid)/population, "chd.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - stroke.preval)/population, "chd.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.cvd.incid - cumul.chd.incid, "cumul.stroke.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - chd.incid)/population, "stroke.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - chd.preval)/population, "stroke.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.chd.deaths, "cumul.chd.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*chd.deaths/population, "chd.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.stroke.deaths, "cumul.stroke.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*stroke.deaths/population, "stroke.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.qalys.disc/1e6, "cumul.qalys.disc.millions", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.disease.cost.disc/1e6, "cumul.disease.cost.disc.millions", prob = c(0.5, 0.025, 0.975))
),
keyby = c("year", "qimd")]
tt <- tt_qimd[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year", "qimd")), qimd+variable~year), output_dir("Report/baseline_summary_results_qimd.csv")) 

# Scenarios
results <- fread(output_dir("results.csv"))
results[, `:=` (
  invitation.cost.disc = invitation.cost * discount ^ (year - 2016),
  participation.cost.disc = participation.cost * discount ^ (year - 2016)
),
by = c("scenario", "mc")]

setkeyv(results, c("year"))
results[, `:=` (
  cumul.invitation.cost.disc = cumsum(invitation.cost.disc),
  cumul.participation.cost.disc = cumsum(participation.cost.disc)
),
by = c("scenario", "mc")]

tt <- results[, c(
  output_quants_lbl(cumul.cpp.cvd, "cumul.cpp.cvd", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.dpp.chd + cumul.dpp.stroke, "cumul.dpp.cvd", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.dpp.other, "cumul.dpp.other", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.net.utility.disc, "cumul.net.utility.disc", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.invitation.cost.disc, "cumul.invitation.cost.disc", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.participation.cost.disc, "cumul.participation.cost.disc", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.net.cost.disc, "cumul.net.cost.disc", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.icer.disc, "cumul.icer.disc", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.nmb.disc, "cumul.nmb.disc", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(sei, "sei", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e6 * rei, "rei", 1, prob = c(0.5, 0.025, 0.975))
),
keyby = c("year", "scenario")]
tt <- tt[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year", "scenario")), scenario+variable~year), output_dir("Report/scenarios_for_report.csv")) 

# graph
tt <- rbind(
  results[, .(type = "Cost-effective", prob = prop_if(cumul.nmb.disc > 0)), by = .(scenario, year)],
  results[, .(type = "Cost saving", prob = prop_if(cumul.net.cost.disc<=0)), by = .(scenario, year)],
  # results[, .(type = "Generate inequalities", prob = prop_if((sei < 0) | (sei > 0 & rei < 0))), by = .(scenario, year)],
  results[, .(type = "Reduce inequalities", prob = prop_if((rei > 0) & (sei > 0))), by = .(scenario, year)]
)
tt[, scenario := factor(scenario)]
tt[, type := factor(type)]
tt[, prob.sm := fbound(predict(loess(prob~year, span = 0.35)), 0, 1), by = .(scenario, type)]


getPalette <- brewer.pal(9, "Set1") 

for (study.sc in levels(tt$scenario)) {
gg <- 
  ggplot(tt[scenario %in% study.sc, ],
         aes(x = year, y = prob.sm, col = type)) +
  geom_hline(aes(yintercept = 0.8),linetype = "dashed", show.legend = F) +
  geom_point(size = 2, alpha = 1/5, show.legend = F) +
  geom_line(size = 1, alpha = 5/5) +
  scale_x_continuous(name="Year", breaks = c(2011,2015, 2020,2025, 2030,2035, 2040)) +
  scale_y_continuous(name="Probability", labels = percent, breaks = seq(0, 1, .2), limits = c(0, 1)) +
  scale_colour_manual(values = getPalette, drop=TRUE, breaks = tt$type, limits = levels(tt$type)) +
  my_theme +
  theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0),
        axis.text.y = element_text(size = rel(0.8), vjust = 0.3, angle = 0))

ggsave(paste0(study.sc, ".tiff"), gg, path = output_dir("Report/"), dpi = 600, scale = 0.7, width = 8, height = 5,
       compression = "lzw", antialias = "subpixel", type="cairo")
save_tikz(paste0(study.sc, ".tex"), gg, path = output_dir("Report/"))
}



# dominant scenarios
results_lim <- results[scenario %in% c("sc00", "sc05a", "sc07b", "sc12b", "sc19", "sc20"), ]
results_lim[scenario == "sc00", scenario := "Current"]
results_lim[scenario == "sc05a", scenario := "Targeted"]
results_lim[scenario == "sc07b", scenario := "Optimal"]
results_lim[scenario == "sc12b", scenario := "Current + Structural"]
results_lim[scenario == "sc19", scenario := "Targeted + Structural"]
results_lim[scenario == "sc20", scenario := "Optimal + Structural"]

nam <- c("Current", "Targeted", "Optimal", "Current + Structural", "Targeted + Structural", "Optimal + Structural")

tt <- rbind(
  results_lim[, .(type = "Cost-effective", prob = prop_if(cumul.nmb.disc > 0)), by = .(scenario, year)],
  results_lim[, .(type = "Cost saving", prob = prop_if(cumul.net.cost.disc<=0)), by = .(scenario, year)],
  # results_lim[, .(type = "Generate inequalities", prob = prop_if((sei < 0) | (sei > 0 & rei < 0))), by = .(scenario, year)],
  results_lim[, .(type = "Reduce inequalities", prob = prop_if((rei > 0) & (sei > 0))), by = .(scenario, year)]
)
tt[, scenario := factor(scenario)]
tt[, type := factor(type)]
tt[, prob.sm := fbound(predict(loess(prob~year, span = 0.35)), 0, 1), by = .(scenario, type)]

getPalette <- brewer.pal(8, "Set2") 

for (type.sc in levels(tt$type)) {
  gg <- 
    ggplot(tt[type %in% type.sc, ],
           aes(x = year, y = prob.sm, col = scenario)) +
    geom_hline(aes(yintercept = 0.8),linetype = "dashed", show.legend = F) +
    geom_point(size = 2, alpha = 1/5, show.legend = F) +
    geom_line(size = 1, alpha = 5/5) +
    scale_x_continuous(name="Year", breaks = c(2011,2015, 2020,2025, 2030,2035, 2040)) +
    scale_y_continuous(name="Probability", labels = percent, breaks = seq(0, 1, .2), limits = c(0, 1)) +
    scale_colour_manual(values = getPalette, drop=TRUE, breaks = nam, limits = levels(tt$scenario)) +
    my_theme +
    theme(axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0),
          axis.text.y = element_text(size = rel(0.8), vjust = 0.3, angle = 0))
  
  ggsave(paste0(type.sc, ".tiff"), gg, path = output_dir("Report/"), dpi = 600, scale = 0.7, width = 8, height = 5,
         compression = "lzw", antialias = "subpixel", type="cairo")
  save_tikz(paste0(type.sc, ".tex"), gg, path = output_dir("Report/"))
}

setkey(tt, year)
fwrite(tt[prob>0.8, utils::head(.SD, 1) , by=.(type, scenario)], output_dir("Report/dominant.csv"))


# current implementation detailed results
validat <- fread(output_dir("results.csv"))
validat <- validat[scenario == "sc00"]
discount <- (100 - 3.5)/100

validat[, `:=` (
  invitation.cost.disc = invitation.cost * discount ^ (year - 2016),
  participation.cost.disc = participation.cost * discount ^ (year - 2016),
  disease.cost.disc = disease.cost * discount ^ (year - 2016),
  qalys.disc = qalys * discount ^ (year - 2016)
)]

setkeyv(validat, c("year"))
validat[, `:=` (
  cumul.chd.incid = cumsum(chd.incid),
  cumul.stroke.incid = cumsum(stroke.incid),
  cumul.cvd.incid = cumsum(cvd.incid),
  cumul.chd.deaths = cumsum(chd.deaths),
  cumul.stroke.deaths = cumsum(stroke.deaths),
  cumul.qalys.disc = cumsum(qalys.disc),
  cumul.disease.cost.disc = cumsum(disease.cost.disc),
  cumul.invitation.cost.disc = cumsum(invitation.cost.disc),
  cumul.participation.cost.disc = cumsum(participation.cost.disc)
),
by = c("scenario", "mc")]


tt <- validat[, c(
  output_quants_lbl(population, "population", prob = c(0.5, 0.025, 0.975), digits = 3),
  output_quants_lbl(cumul.cvd.incid - cumul.stroke.incid, "cumul.chd.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - stroke.incid)/population, "chd.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - stroke.preval)/population, "chd.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.cvd.incid - cumul.chd.incid, "cumul.stroke.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - chd.incid)/population, "stroke.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - chd.preval)/population, "stroke.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.chd.deaths, "cumul.chd.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*chd.deaths/population, "chd.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.stroke.deaths, "cumul.stroke.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*stroke.deaths/population, "stroke.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.qalys.disc/1e6, "cumul.qalys.disc.millions", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.disease.cost.disc/1e6, "cumul.disease.cost.disc.millions", prob = c(0.5, 0.025, 0.975))
),
keyby = c("year")]
tt <- tt[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year")), variable~year), output_dir("Report/current_summary_results.csv")) 



validat <- fread(output_dir("results_sex.csv"))
validat <- validat[scenario == "sc00"]
discount <- (100 - 3.5)/100

validat[, `:=` (
  invitation.cost.disc = invitation.cost * discount ^ (year - 2016),
  participation.cost.disc = participation.cost * discount ^ (year - 2016),
  disease.cost.disc = disease.cost * discount ^ (year - 2016),
  qalys.disc = qalys * discount ^ (year - 2016)
)]

setkeyv(validat, c("year"))
validat[, `:=` (
  cumul.chd.incid = cumsum(chd.incid),
  cumul.stroke.incid = cumsum(stroke.incid),
  cumul.cvd.incid = cumsum(cvd.incid),
  cumul.chd.deaths = cumsum(chd.deaths),
  cumul.stroke.deaths = cumsum(stroke.deaths),
  cumul.qalys.disc = cumsum(qalys.disc),
  cumul.disease.cost.disc = cumsum(disease.cost.disc),
  cumul.invitation.cost.disc = cumsum(invitation.cost.disc),
  cumul.participation.cost.disc = cumsum(participation.cost.disc)
),
by = c("scenario", "mc", "sex")]

tt_sex <- validat[, c(
  output_quants_lbl(population, "population", prob = c(0.5, 0.025, 0.975), digits = 3),
  output_quants_lbl(cumul.cvd.incid - cumul.stroke.incid, "cumul.chd.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - stroke.incid)/population, "chd.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - stroke.preval)/population, "chd.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.cvd.incid - cumul.chd.incid, "cumul.stroke.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - chd.incid)/population, "stroke.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - chd.preval)/population, "stroke.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.chd.deaths, "cumul.chd.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*chd.deaths/population, "chd.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.stroke.deaths, "cumul.stroke.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*stroke.deaths/population, "stroke.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.qalys.disc/1e6, "cumul.qalys.disc.millions", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.disease.cost.disc/1e6, "cumul.disease.cost.disc.millions", prob = c(0.5, 0.025, 0.975))
),
keyby = c("year", "sex")]
tt <- tt_sex[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year", "sex")), sex+variable~year), output_dir("Report/current_summary_results_sex.csv")) 


validat <- fread(output_dir("results_age.csv"))
validat <- validat[scenario == "sc00"]
discount <- (100 - 3.5)/100

validat[, `:=` (
  invitation.cost.disc = invitation.cost * discount ^ (year - 2016),
  participation.cost.disc = participation.cost * discount ^ (year - 2016),
  disease.cost.disc = disease.cost * discount ^ (year - 2016),
  qalys.disc = qalys * discount ^ (year - 2016)
)]

setkeyv(validat, c("year"))
validat[, `:=` (
  cumul.chd.incid = cumsum(chd.incid),
  cumul.stroke.incid = cumsum(stroke.incid),
  cumul.cvd.incid = cumsum(cvd.incid),
  cumul.chd.deaths = cumsum(chd.deaths),
  cumul.stroke.deaths = cumsum(stroke.deaths),
  cumul.qalys.disc = cumsum(qalys.disc),
  cumul.disease.cost.disc = cumsum(disease.cost.disc),
  cumul.invitation.cost.disc = cumsum(invitation.cost.disc),
  cumul.participation.cost.disc = cumsum(participation.cost.disc)
),
by = c("scenario", "mc", "agegroup")]

tt_age <- validat[, c(
  output_quants_lbl(population, "population", prob = c(0.5, 0.025, 0.975), digits = 3),
  output_quants_lbl(cumul.cvd.incid - cumul.stroke.incid, "cumul.chd.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - stroke.incid)/population, "chd.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - stroke.preval)/population, "chd.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.cvd.incid - cumul.chd.incid, "cumul.stroke.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - chd.incid)/population, "stroke.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - chd.preval)/population, "stroke.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.chd.deaths, "cumul.chd.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*chd.deaths/population, "chd.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.stroke.deaths, "cumul.stroke.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*stroke.deaths/population, "stroke.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.qalys.disc/1e6, "cumul.qalys.disc.millions", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.disease.cost.disc/1e6, "cumul.disease.cost.disc.millions", prob = c(0.5, 0.025, 0.975))
),
keyby = c("year", "agegroup")]
tt <- tt_age[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year", "agegroup")), agegroup+variable~year), output_dir("Report/current_summary_results_age.csv")) 



validat <- fread(output_dir("results_qimd.csv"))
validat <- validat[scenario == "sc00"]
discount <- (100 - 3.5)/100

validat[, `:=` (
  invitation.cost.disc = invitation.cost * discount ^ (year - 2016),
  participation.cost.disc = participation.cost * discount ^ (year - 2016),
  disease.cost.disc = disease.cost * discount ^ (year - 2016),
  qalys.disc = qalys * discount ^ (year - 2016)
)]

setkeyv(validat, c("year"))
validat[, `:=` (
  cumul.chd.incid = cumsum(chd.incid),
  cumul.stroke.incid = cumsum(stroke.incid),
  cumul.cvd.incid = cumsum(cvd.incid),
  cumul.chd.deaths = cumsum(chd.deaths),
  cumul.stroke.deaths = cumsum(stroke.deaths),
  cumul.qalys.disc = cumsum(qalys.disc),
  cumul.disease.cost.disc = cumsum(disease.cost.disc),
  cumul.invitation.cost.disc = cumsum(invitation.cost.disc),
  cumul.participation.cost.disc = cumsum(participation.cost.disc)
),
by = c("scenario", "mc", "qimd")]

tt_qimd <- validat[, c(
  output_quants_lbl(population, "population", prob = c(0.5, 0.025, 0.975), digits = 3),
  output_quants_lbl(cumul.cvd.incid - cumul.stroke.incid, "cumul.chd.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - stroke.incid)/population, "chd.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - stroke.preval)/population, "chd.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.cvd.incid - cumul.chd.incid, "cumul.stroke.incid", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.incid - chd.incid)/population, "stroke.incid.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5 * (cvd.preval - chd.preval)/population, "stroke.preval.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.chd.deaths, "cumul.chd.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*chd.deaths/population, "chd.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.stroke.deaths, "cumul.stroke.deaths", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(1e5*stroke.deaths/population, "stroke.deaths.rate", 1, prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.qalys.disc/1e6, "cumul.qalys.disc.millions", prob = c(0.5, 0.025, 0.975)),
  output_quants_lbl(cumul.disease.cost.disc/1e6, "cumul.disease.cost.disc.millions", prob = c(0.5, 0.025, 0.975))
),
keyby = c("year", "qimd")]
tt <- tt_qimd[year %in% c(2011, 2020, 2030, 2040), ]
fwrite(dcast(melt(tt, id.vars = c("year", "qimd")), qimd+variable~year), output_dir("Report/current_summary_results_qimd.csv")) 
