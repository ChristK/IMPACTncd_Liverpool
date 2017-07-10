#cmpfile("./post_simulation_functions.R")
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
dependencies(c("Kmisc", 
               "ggplot2",
               "scales",
               "compiler",
               "survey",
               "ggrepel",
               "doParallel",
               "doRNG",
               "foreach",
               "future",
               "data.table"))



# get object name. Doesn't work inside a local function when a set of list items are passed to lapply (and it also fails when an object is passed from a list given to a for-loop.) 
get_object_name <- 
    function(v1) { 
      s <- substitute(v1)
      if (length(s) == 1) deparse(s) else sub("\\(.", "", s[2])
    }

# output quantiles
output_quants <- 
  cmpfun(
    function(x, nam = get_object_name(x), pop_extrapol = pop.fraction, ...) {
      pr <- c(0.025, .2, 0.25, 0.5, 0.75, 0.8, 0.975)
      pr_nam <- percent(pr) #Very Slow
      out <- transpose(list(round(quantile(x, pr, T,
                                    F) / pop_extrapol, ...)))
      setattr(out, "names", paste0(nam, "_", pr_nam))
      return(out)
    }
  )

qimd_labeller <- labeller(qimd = c("1" = "QIMD 1\nleast deprived",
                                   "2" = "QIMD 2",
                                   "3" = "QIMD 3",
                                   "4" = "QIMD 4",
                                   "5" = "QIMD 5\nmost deprived"))
gbp <- dollar_format(prefix = "", suffix = "£")

process_results <- 
  cmpfun(
    function(dt, discount_pct = 3.5, strat = NULL) {
      discount <- (100 - discount_pct)/100
     
      dt[, `:=` (
        qalys.disc = qalys * discount ^ (year - 2016),
        net.utility.disc = net.utility * discount ^ (year - 2016),
        net.utility.curr.impl.disc = net.utility.curr.impl * discount ^ (year - 2016),
        net.cost.disc = net.cost * discount ^ (year - 2016),
        net.cost.curr.impl.disc = net.cost.curr.impl * discount ^ (year - 2016)
      ),
      by = c("scenario", "mc", strat)]
      
      setkeyv(dt, c("year", strat))
      dt[, `:=` (
        cumul.qalys.disc = cumsum(qalys.disc),
        cumul.net.utility.disc = cumsum(net.utility.disc),
        cumul.net.utility.curr.impl.disc = cumsum(net.utility.curr.impl.disc),
        cumul.net.cost.disc = cumsum(net.cost.disc),
        cumul.net.cost.curr.impl.disc = cumsum(net.cost.curr.impl.disc),
        cumul.net.utility = cumsum(net.utility),
        cumul.net.utility.curr.impl = cumsum(net.utility.curr.impl),
        cumul.net.cost = cumsum(net.cost),
        cumul.net.cost.curr.impl = cumsum(net.cost.curr.impl),
        cumul.cpp.chd = cumsum(cpp.chd),
        cumul.cpp.chd.curr.impl = cumsum(cpp.chd.curr.impl),
        cumul.cpp.stroke = cumsum(cpp.stroke),
        cumul.cpp.stroke.curr.impl = cumsum(cpp.stroke.curr.impl),
        cumul.cpp.cvd = cumsum(cpp.cvd),
        cumul.cpp.cvd.curr.impl = cumsum(cpp.cvd.curr.impl),
        cumul.dpp.chd = cumsum(dpp.chd),
        cumul.dpp.chd.curr.impl = cumsum(dpp.chd.curr.impl),
        cumul.dpp.stroke = cumsum(dpp.stroke),
        cumul.dpp.stroke.curr.impl = cumsum(dpp.stroke.curr.impl),
        cumul.dpp.other = cumsum(dpp.other),
        cumul.dpp.other.curr.impl = cumsum(dpp.other.curr.impl)
      ),
      by = c("scenario", "mc", strat)]
      
      # dt[net.utility.disc == 0, net.utility.disc := 0.1]
      # dt[net.utility == 0, net.utility := 0.1]
      # dt[cumul.net.utility.disc == 0, cumul.net.utility.disc := 0.1]
      # dt[cumul.net.utility == 0, cumul.net.utility := 0.1]
      # dt[net.utility.curr.impl.disc == 0, net.utility.curr.impl.disc := 0.1]
      # dt[net.utility.curr.impl == 0, net.utility.curr.impl := 0.1]
      # dt[cumul.net.utility.curr.impl.disc == 0, cumul.net.utility.curr.impl.disc := 0.1]
      # dt[cumul.net.utility.curr.impl == 0, cumul.net.utility.curr.impl := 0.1]

      dt[
        ,`:=`(icer.disc       = net.cost.disc / net.utility.disc,
              nmb.disc        = -net.cost.disc + net.utility.disc * 20000,
              cumul.icer.disc = cumul.net.cost.disc / cumul.net.utility.disc,
              cumul.nmb.disc  = -cumul.net.cost.disc + cumul.net.utility.disc * 20000,
              icer.curr.impl.disc       = net.cost.curr.impl.disc / net.utility.curr.impl.disc,
              nmb.curr.impl.disc        = -net.cost.curr.impl.disc + net.utility.curr.impl.disc * 20000,
              cumul.icer.curr.impl.disc = cumul.net.cost.curr.impl.disc / cumul.net.utility.curr.impl.disc,
              cumul.nmb.curr.impl.disc  = -cumul.net.cost.curr.impl.disc + cumul.net.utility.curr.impl.disc * 20000,
              icer       = net.cost / net.utility,
              nmb        = -net.cost + net.utility * 20000,
              cumul.icer = cumul.net.cost / cumul.net.utility,
              cumul.nmb  = -cumul.net.cost + cumul.net.utility * 20000,
              icer.curr.impl       = net.cost.curr.impl / net.utility.curr.impl,
              nmb.curr.impl        = -net.cost.curr.impl + net.utility.curr.impl * 20000,
              cumul.icer.curr.impl = cumul.net.cost.curr.impl / cumul.net.utility.curr.impl,
              cumul.nmb.curr.impl  = -cumul.net.cost.curr.impl + cumul.net.utility.curr.impl * 20000
        )]
      # replace non finite
      if (!is.null(strat)) {strat <- ifelse(strat == "age", "agegroup", strat)}
      
      for (j in outersect(names(dt), c(strat, "year", "scenario", "mc")))
        set(dt,which(!is.finite(dt[[j]])),j,NA)
    }
  )

summarise_results <-
  cmpfun(function(dt,
                  probs = c(0.5, 0.025, 0.975, .2, 0.8, 0.25, 0.75),
                  pop_extrapol = pop.fraction,
                  digits = 0,
                  strat = NULL,
                  ...) {
    # define variables that need not to be divided by pop.fraction
    protect <- c("icer", ".prevalence", ".5ormore", ".mean")
    out <- suppressWarnings(melt(dt, id.vars = c("year", "mc", "scenario", strat)))
    out[is.infinite(value), value := NA]
    out[, 
        id:=.GRP,
        keyby = c("year", strat, "scenario", "variable")]
    setkey(out, id) # Need to be before fquantile_byid
    out[, paste0("V", seq_along(probs)) := mat2df(fquantile_byid(value, probs, id))]
    out2 <- out[, setDT(mat2df(fquantile_byid(value, probs, id)))]
    out2[, id := .I] # because out was sorted on id
    out <- unique(out[, .SD, .SDcols = c("id", "year", strat, "scenario", "variable")])
    out <- merge(out, out2, by = "id", all = T)
    out[, id := NULL]
    setkeyv(out, c("scenario", "variable", strat, "year"))
    
    # select lines that variable name does not include any of the protected terms
    # and divide by pop_extrapol
    out[Reduce(`+`, lapply(protect, grepl, variable)) == 0, 
        paste0("V", seq_along(probs)) := 
          lapply(mget(paste0("V", seq_along(probs))),
                 function(x) x/pop_extrapol)]
    # apply rounding
    out[,
        paste0("V", seq_along(probs)) :=
          lapply(mget(paste0("V", seq_along(probs))), round, digits)]
    setnames(out, paste0("V", seq_along(probs)), percent(probs))
    
    # remove lines with NA
    out <- na.omit(out)
    return(out)
  })

my_theme <-   theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = rel(1), vjust = 0),
        axis.title = element_text(size = rel(1)),
        axis.title.x = element_text(margin = margin(20,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(size = rel(0.8), vjust = 0, angle = 0),
        axis.text.y = element_text(size = rel(0.8), vjust = 0),
        strip.text.x = element_text(rel(1.2)),
        strip.text.y = element_text(rel(1.2)),
        strip.background = element_rect(colour = "lavenderblush4", fill = "lavenderblush2"),
        legend.title = element_blank(), # no legent title
        legend.position = "right",
        #legend.justification = c(0,1), legend.position = c(0,1),
        legend.key = element_rect(size = 1, colour = "transparent", fill = "transparent"),
        legend.key.size = unit(1, "lines")) 

# Remove iteration files that their log is not the same length as the most frequent log
remove_faulty_iterations <-
  function(path_to_output = output_dir("Intermediate"))
  {
    all.files <- as.list(dir(
      path = path_to_output,
      pattern = glob2rx("*_log.txt"),
      full.names = T,
      recursive = F
    ))
    if (length(all.files) > 0) {
    logs <- mclapply(all.files,
                     readLines,
                     mc.cores = clusternumber)
    length(logs)
    lapply(logs, length)
    # get the most common length. max length doesn't work because crashes sometimes fill log with gibberish
    common_length <- which.max(tabulate(sapply(logs, length)))
    
    common_length <- which(sapply(logs, function(x) {
      ifelse(length(x) == common_length, FALSE, TRUE)
    }))
    all.files <-
      unlist(all.files[common_length]) # only logs of non common length
    
    all.files <- gsub(path_to_output, "", all.files, fixed = T)
    
    all.files <- gsub("_log.txt", "", all.files, fixed = T)
    
    all.files2 <- lapply(all.files, function(x) {
      dir(
        path = path_to_output,
        pattern = x,
        full.names = T,
        recursive = F
      )
    })
    
    lapply(all.files2, file.remove)
    }
    return(NULL)
  }

# Ensures that when fwrites appent file colnames of file to be written, match those already in the file
fwrite_safe <- function(x,
                        file = "",
                        append = FALSE,
                        ...)
{
  if (append == TRUE)
  {
    if (file.exists(file)) {
      col_names_disk <- names(fread(file, nrows = 0))
      col_names_file <- names(x)
      col_names <- outersect(col_names_disk, col_names_file)
      if (length(col_names) > 0)
        x[, (col_names) := NA]
      setcolorder(x, col_names_disk)
    }
  }
  fwrite(x, file, append, ...)
}

