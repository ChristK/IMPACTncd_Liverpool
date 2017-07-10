#cmpfile("./cluster functions.R")
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

haha <- paste0(sample(c(rep(0:9, each = 5), LETTERS, letters), 12,
                      replace = T),
               collapse = "")

# Define function for output dir
output.dir <-
  cmpfun(function(x = NULL) {
    output_dir(paste0("Intermediate/", haha, "_", x))
  })

# OLD function
# output.dir <-
#   cmpfun(function(x = NULL) {
#     paste0("./Output/Intermediate/", haha, "/", x)
#   })

dir.create(path = output_dir("Intermediate/"),
           recursive = T,
           showWarnings = F) # create a unique directory for each run of each scenario

# from Mozaffarian NEJM
salt.sbp.reduct <-
  cmpfun(function(salt.difference, age, sbp, race, n) {
    y = rnorm(n, -3.735, 0.73) + rnorm(n, -0.105, 0.029) * (age - 50) +
      rnorm(n, -1.874, 0.884) * (sbp > 140) +
      rnorm(n, -2.489, 1.188) * (race == 6)
    return(salt.difference * y / 5.85)
  })

# Define function to match continuous distributions of each group with the one in SPOP2011 to simulate ageing
# ageing.distr <- cmpfun(function(risk.factor, env = my.env) {
#   risk.factor <- substitute(risk.factor)
#   temp <- SPOP2011[eval(risk.factor)>0, list(eval(risk.factor), group)]
#   nam <- paste0(risk.factor, ".rank")
#   temp[, (nam) := (frank(eval(risk.factor), na.last = F, ties.method="random")-1)/(.N - 1),
#        by = group]
#   setkeyv(temp, c("group", nam))
#
#   POP[, (nam) := (frank(eval(risk.factor), na.last = F, ties.method="random")-1)/(.N - 1),
#       by = group]
#   POP[, eval(risk.factor) := NULL]
#   setkeyv(POP, c("group", nam))
#   POP <- temp[POP, roll = "nearest"]
#   assign("POP", POP, envir = env)
# }
# )

ageing.distr <-
  # smaller fortune increase the variability of the join
  cmpfun(function(DT_ref, risk_factor, fortune = 50) {
    nam  <- paste0(substitute(risk_factor), ".rank")
    POP[, (risk_factor) := NULL]
    setkeyv(POP, c("group", nam))
    min_group <- min(DT_ref[, .N, by = group][, min(N)], fortune)
    tt <-
      DT_ref[, .SD[sample(.N, min_group)], by = group] # +min_group otherwise small groups get 0
    setkeyv(tt, c("group", nam))
    assign("POP", tt[POP, roll = "nearest"], envir = my.env)
  })

prepare_outcomes <-
  cmpfun(function(dt) {
    outcome_names <-
      c(
        "chd.incidence",
        "stroke.incidence",
        "other.mortality",
        "chd.mortality",
        "stroke.mortality",
        "htn.incidence",
        "diabtotr"
      )
    
    # Logic to only collect scenario RF if present in the DT
    sc_names <-   gsub(
      ".R",
      "",
      list.files(
        path = "./Scenarios",
        pattern = glob2rx("sc*.R"),
        full.names = F,
        recursive = F
      )
    )
    sc_outcome_names <-
      unlist(lapply(outcome_names, function(outcome_names, sc_names)
        paste(outcome_names, sc_names, sep = "_"), sc_names),
        use.names = F)
    sc_outcome_names <- intersect(sc_outcome_names, names(dt))
    
    if (length(sc_outcome_names) == 0) {
      nam <- intersect(outcome_names, names(dt))
    } else
      nam <- sc_outcome_names
    
    if (!all(key(dt) == c("id", "year")))
      setkey(dt, id, year)
    dt[, diabtotr_out := as.integer(get(grep("diabtotr", nam, value = T))) - 1L]
    dt[diabtotr_out == 1L, diabtotr_out := 700L]
    dt[, chd.incidence_out := get(grep("chd.incidence", nam, value = T))]
    dt[chd.incidence_out > 1L, chd.incidence_out := 400000L]
    dt[chd.incidence_out == 1L, chd.incidence_out := 400001L] # Only 1st year important for costs
    dt[, stroke.incidence_out := get(grep("stroke.incidence", nam, value = T))]
    dt[stroke.incidence_out > 1L, stroke.incidence_out := 50000L]
    dt[stroke.incidence_out == 1L, stroke.incidence_out := 50010L] # Only 1st year important for costs
    
    t1 <-
      df2mat(dcast(dt[, .(id, year, get(grep("other.mortality", nam, value = T)))][year >= 0, ],
                   id ~ year, value.var = "V3")) # -100
    t2 <-
      df2mat(dcast(dt[, .(id, year, get(grep("chd.mortality", nam, value = T)))][year >= 0, ],
                   id ~ year, value.var = "V3")) # -20
    t3 <-
      df2mat(dcast(dt[, .(id, year, get(grep("stroke.mortality", nam, value = T)))][year >= 0, ],
                   id ~ year, value.var = "V3")) # -3
    t4 <-
      df2mat(dcast(dt[, .(id, year, chd.incidence_out)][year >= 0, ],
                   id ~ year, value.var = "chd.incidence_out")) # 400000
    t5 <-
      df2mat(dcast(dt[, .(id, year, stroke.incidence_out)][year >= 0, ],
                   id ~ year, value.var = "stroke.incidence_out")) # 50000
    t6 <-
      df2mat(dcast(dt[, .(id, year, get(grep("htn.incidence", nam, value = T)))][year >= 0, ],
                   id ~ year, value.var = "V3")) # 6000
    t7 <- df2mat(dcast(dt[, .(id, year, diabtotr_out)][year >= 0, ],
                       id ~ year, value.var = "diabtotr_out")) # 700
    
    a0 <- df2mat(dcast(dt[, .(id, year, age)][year >= 0, ],
                       id ~ year, value.var = "age"))
    a1 <-
      df2mat(utility.mc[condition == "event_free", .(age, qaly)])
    a2 <- matrix(
      c(
        400000,
        50000,
        6000,
        700,
        utility.mc[condition == "chd_exclusive", unique(qaly)],
        utility.mc[condition == "stroke_exclusive", unique(qaly)],
        utility.mc[condition == "htn_exclusive", unique(qaly)],
        utility.mc[condition == "diabetes_exclusive", unique(qaly)]
      ),
      nrow = 4,
      ncol = 2
    )
    a3 <-
      df2mat(dcast(dt[, .(id, year, qimd)][year >= 0, ], id ~ year, value.var = "qimd"))
    a4 <-
      dcast(cost.mc, year + qimd ~ condition, value.var = "cost")[, `:=` (chd_death = chd.cost.death.mc,
                                                                          stroke_death = stroke.cost.death.mc,
                                                                          qimd = as.integer(qimd))]
    a4 <- df2mat(a4)
     out <- collect_output(t1, t2, t3, t4, t5, t6, t7, a0, a1, a2, a3, a4)
     
     outcome_names <-
       c(
         "invited",
         "participated"
       )
     sc_outcome_names <-
       unlist(lapply(outcome_names, function(outcome_names, sc_names)
         paste(outcome_names, sc_names, sep = "_"), sc_names),
         use.names = F)
     sc_outcome_names <- intersect(sc_outcome_names, names(dt))
     
     if (length(sc_outcome_names) > 0) {
       out$invited <- df2mat(dcast(dt[, .(id, year, get(grep("invited", sc_outcome_names, value = T)))][year>=0,],
                                   id~year, value.var = "V3")) 
       
       out$participated <- df2mat(dcast(dt[, .(id, year, get(grep("participated", sc_outcome_names, value = T)))][year>=0,],
                                        id~year, value.var = "V3"))
       
       c3 <- which(is.na(out$lifecourse))
       out$invited[c3] <- NA
       out$participated[c3] <- NA
     }
    return(out)
  })

summary_fn <-
  cmpfun(function(dt, nam, sc_name, strat = NULL) {
    output <- cbind.dt(
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          count_if(x >= 0)),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "population"
      ),
      melt(
        dt$lifecourse_masked65[, lapply(.SD, function(x)
          count_if(x >= 0)),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "population30to65"
      )[, .(population30to65)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          # exclusively. ie chd + stroke is not counted here
          count_if((x %% 100L) %in% c(1))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "chd.incid"
      )[, .(chd.incid)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          # exclusively. ie chd + stroke is not counted here
          count_if((x %% 100L) %in% c(10))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "stroke.incid"
      )[, .(stroke.incid)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          count_if((x %% 100L) %in% c(1, 10, 11))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "cvd.incid"
      )[, .(cvd.incid)],
      
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          # exclusively. ie chd + stroke is not counted here
          count_if(between(x, 400000, 410000))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "chd.preval"
      )[, .(chd.preval)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          # exclusively. ie chd + stroke is not counted here
          count_if(between(x, 50000, 60000))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "stroke.preval"
      )[, .(stroke.preval)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          count_if(between(x, 50000, 500000))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "cvd.preval"
      )[, .(cvd.preval)],
      
      melt(
        dt$lifecourse_masked65[, lapply(.SD, function(x)
          # exclusively. ie chd + stroke is not counted here
          count_if(between(x, 400000, 410000))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "chd.preval.under65"
      )[, .(chd.preval.under65)],
      melt(
        dt$lifecourse_masked65[, lapply(.SD, function(x)
          # exclusively. ie chd + stroke is not counted here
          count_if(between(x, 50000, 60000))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "stroke.preval.under65"
      )[, .(stroke.preval.under65)],
      melt(
        dt$lifecourse_masked65[, lapply(.SD, function(x)
          count_if(between(x, 50000, 500000))),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "cvd.preval.under65"
      )[, .(cvd.preval.under65)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          count_if(x == -20)),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "chd.deaths"
      )[, .(chd.deaths)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          count_if(x == -3)),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "stroke.deaths"
      )[, .(stroke.deaths)],
      melt(
        dt$lifecourse[, lapply(.SD, function(x)
          count_if(x == -100)),
          .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "other.deaths"
      )[, .(other.deaths)],
      melt(
        dt$cost[, lapply(.SD, sum, na.rm = T), .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "cost"
      )[, .(disease.cost = cost)],
      switch(("invited" %in% names(dt)) + 1L,
             NULL,
             melt(
               dt$invited[, lapply(.SD, sum, na.rm = T),
                          .SDcols = nam, keyby = strat],
               measure.vars = nam,
               variable.name = "year",
               value.name = "cost"
             )[, .(invitation.cost = invitation_cost * cost)]
      ),
      switch(("participated" %in% names(dt)) + 1L,
             NULL,
             melt(
               dt$participated[, lapply(.SD, sum, na.rm = T),
                               .SDcols = nam, keyby = strat],
               measure.vars = nam,
               variable.name = "year",
               value.name = "cost"
             )[, .(participation.cost = participation_cost * cost)]
      ),
      melt(
        dt$utility[, lapply(.SD, sum, na.rm = T), .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "qalys"
      )[, .(qalys)]
    )
     if (is.null(strat)) {
    if (exists("out_summary")) {
      output[out_summary,
             `:=`(net.utility = qalys - i.qalys,
                  net.cost    = disease.cost + participation.cost + invitation.cost - i.disease.cost,
                  cpp.chd     = i.chd.incid - chd.incid,
                  cpp.stroke  = i.stroke.incid - stroke.incid,
                  cpp.cvd     = i.cvd.incid - cvd.incid,
                  dpp.chd     = i.chd.deaths - chd.deaths,
                  dpp.stroke  = i.stroke.deaths - stroke.deaths,
                  dpp.other   = i.other.deaths - other.deaths),
             on = "year"]
    }
       if (exists("out_summary_sc00")) {
         output[
           out_summary_sc00,
           `:=`(net.utility.curr.impl    = qalys - i.qalys,
                net.cost.curr.impl       = disease.cost + participation.cost + invitation.cost -
                  i.disease.cost - i.participation.cost - i.invitation.cost,
                cpp.chd.curr.impl    = i.chd.incid - chd.incid,
                cpp.stroke.curr.impl = i.stroke.incid - stroke.incid,
                cpp.cvd.curr.impl    = i.cvd.incid - cvd.incid,
                dpp.chd.curr.impl    = i.chd.deaths - chd.deaths,
                dpp.stroke.curr.impl = i.stroke.deaths - stroke.deaths,
                dpp.other.curr.impl  = i.other.deaths - other.deaths),
           on = "year"]
       } else {
         if (exists("out_summary")) {
           output[out_summary,
                  `:=`(net.utility.curr.impl      = 0,
                       net.cost.curr.impl         = 0,
                       cpp.chd.curr.impl          = 0,
                       cpp.stroke.curr.impl       = 0,
                       cpp.cvd.curr.impl          = 0,
                       dpp.chd.curr.impl          = 0,
                       dpp.stroke.curr.impl       = 0,
                       dpp.other.curr.impl        = 0),
                  on = "year"]
         }
       }
     } else { #if strat not NULL
       if (strat == "age") {
         output[, age := age - init.year + as.integer(factor_to_char(year))] # align age with year
         output <- output[between(age, ageL, ageH)]
         output <- agegroup_aggregation(output)
       }
       if (exists(paste0("out_summary", "_", strat))) {
         output[get(paste0("out_summary", "_", strat)),
                `:=`(net.utility = qalys - i.qalys,
                     net.cost    = disease.cost + participation.cost + invitation.cost - i.disease.cost,
                     cpp.chd     = i.chd.incid - chd.incid,
                     cpp.stroke  = i.stroke.incid - stroke.incid,
                     cpp.cvd     = i.cvd.incid - cvd.incid,
                     dpp.chd     = i.chd.deaths - chd.deaths,
                     dpp.stroke  = i.stroke.deaths - stroke.deaths,
                     dpp.other   = i.other.deaths - other.deaths),
                on = c("year", ifelse(strat == "age", "agegroup", strat))]
       }
       if (exists(paste0("out_summary", "_", strat, "_sc00"))) {
         output[
           get(paste0("out_summary", "_", strat, "_sc00")),
           `:=`(net.utility.curr.impl    = qalys - i.qalys,
                net.cost.curr.impl       = disease.cost + participation.cost + invitation.cost -
                  i.disease.cost - i.participation.cost - i.invitation.cost,
                cpp.chd.curr.impl    = i.chd.incid - chd.incid,
                cpp.stroke.curr.impl = i.stroke.incid - stroke.incid,
                cpp.cvd.curr.impl    = i.cvd.incid - cvd.incid,
                dpp.chd.curr.impl    = i.chd.deaths - chd.deaths,
                dpp.stroke.curr.impl = i.stroke.deaths - stroke.deaths,
                dpp.other.curr.impl  = i.other.deaths - other.deaths),
           on = c("year", ifelse(strat == "age", "agegroup", strat))]
       } else {
         if (exists("out_summary")) {
           output[out_summary,
                  `:=`(net.utility.curr.impl      = 0L,
                       net.cost.curr.impl         = 0L,
                       cpp.chd.curr.impl          = 0L,
                       cpp.stroke.curr.impl       = 0L,
                       cpp.cvd.curr.impl          = 0L,
                       dpp.chd.curr.impl          = 0L,
                       dpp.stroke.curr.impl       = 0L,
                       dpp.other.curr.impl        = 0L),
                  on = "year"]
         }
       }
       }
    
    set(output, NULL, "mc", haha)
    set(output, NULL, "scenario", sc_name)
    
    return(output)
  })

mask_age <- cmpfun(function(l) {
  nam <- paste0(init.year:(init.year + yearstoproject - 1))
  
  l <- lapply(l, data.table)
  l <- lapply(l, setnames,
              c("id", nam))
  l <- lapply(l, merge, index, by = "id")
  
  # mask not needed ages
  l$lifecourse_masked65 = copy(l$lifecourse)
  
  for (jj in seq_along(nam))
  {
    lapply(l, function(x) {
      set(x,
          which(!between(index$age , ageL - jj + 1L , ageH - jj + 1L)),
          nam[[jj]], NA)
    })
    set(l$lifecourse_masked65,
        which(!between(index$age , 30 - jj + 1L , 64 - jj + 1L)),
        nam[[jj]], NA)
  }
  return(l)
})

# Collect RF from POP (it has to run before out converted to list of DT)
collect_RF <- cmpfun(function(dt) {
  RF_names <-
    c(
      "cigst1",
      "expsmokCat",
      "a30to06m",
      "porftvg",
      "omsysval",
      "cholval",
      "bmival",
      "diabtotr"
    )
  
  # Logic to only collect scenario RF if present in the DT
  sc_names <-   gsub(
    ".R",
    "",
    list.files(
      path = "./Scenarios",
      pattern = glob2rx("sc*.R"),
      full.names = F,
      recursive = F
    )
  )
  sc_RF_names <-
    unlist(lapply(RF_names, function(RF_names, sc_names)
      paste(RF_names, sc_names, sep = "_"), sc_names), use.names = F)
  sc_RF_names <- intersect(sc_RF_names, names(dt))
  if (length(sc_RF_names) == 0) {
    nam <- c("age", intersect(RF_names, names(dt)))
  } else
    nam <- c("age", sc_RF_names)
  
  dead_position <-
    which(is.na(out$lifecourse)) # need to run before convertion to DT
  output <- vector("list", length(nam))
  names(output) <- nam
  output <- lapply(nam, function(x, output, dead_position) {
    tt <-  df2mat(dcast(dt[, .(id, year, get(x))][year >= 0,],
                        id ~ year, value.var = "V3"))
    tt[dead_position] <-
      NA # Replace risk factor for deads with NAs.
    output[[x]] <- tt
  }, output, dead_position)
  
  setattr(output, "names", nam) # Avoid any shallow copy
  
  ## convert to DT
  namyear <- paste0(init.year:(init.year + yearstoproject - 1))
  output <- lapply(output, data.table)
  output <- lapply(output, setnames,
                   c("id", namyear))
  output <- lapply(output, merge, index, by = "id")
  
  # mask  not needed ages
  for (jj in seq_along(namyear))
  {
    lapply(output, function(x) {
      set(x,
          which(!between(index$age , ageL - jj + 1L , ageH - jj + 1L)),
          namyear[[jj]], NA)
    })
  }
  
  # rename to static names
  setattr(output, "names", gsub("\\_.*", "", names(output)))
  return(output)
})


# Export RF
summary_RF_fn <-
  cmpfun(function(dt, nam, sc_name, strat = NULL) {
    output <- cbind.dt(
      melt(
        dt$age[, lapply(.SD, function(x)
          count_if(x >= 0)), .SDcols = nam, keyby = strat],
        measure.vars = nam,
        variable.name = "year",
        value.name = "population"
      ),
      if ("cigst1" %in% names(dt))
        melt(
          dt$cigst1[, lapply(.SD, function(x)
            prop_if(x == 4)), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "smoking.prevalence"
        ),
      if ("cigst1" %in% names(dt))
        melt(
          dt$cigst1[, lapply(.SD, function(x)
            prop_if(x %in% c(2, 3))), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "exsmoking.prevalence"
        ),
      if ("expsmokCat" %in% names(dt))
        melt(
          dt$expsmokCat[, lapply(.SD, function(x)
            prop_if(x == 2)), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "ets.prevalence"
        ),
      if ("a30to06m" %in% names(dt))
        melt(
          dt$a30to06m[, lapply(.SD, function(x)
            prop_if(x > 4)), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "pa.prevalence.5ormore"
        ),
      if ("porftvg" %in% names(dt))
        melt(
          dt$porftvg[, lapply(.SD, function(x)
            prop_if(x > 4)), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "fv.prevalence.5ormore"
        ),
      if ("diabtotr" %in% names(dt))
        melt(
          dt$diabtotr[, lapply(.SD, function(x)
            prop_if(x == 2)), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "diab.prevalence"
        ),
      if ("bmival" %in% names(dt))
        melt(
          dt$bmival[, lapply(.SD, mean, na.rm = T), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "bmi.mean"
        ),
      if ("omsysval" %in% names(dt))
        melt(
          dt$omsysval[, lapply(.SD, mean, na.rm = T), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "sbp.mean"
        ),
      if ("cholval" %in% names(dt))
        melt(
          dt$cholval[, lapply(.SD, mean, na.rm = T), .SDcols = nam, keyby = strat],
          measure.vars = nam,
          variable.name = "year",
          value.name = "tc.mean"
        )
    )
    output <- output[, .SD, .SDcols = (unique(names(output)))]
    set(output, NULL, "mc", haha)
    set(output, NULL, "scenario", sc_name)
    
    if (!is.null(strat) && strat == "age") {
      output[, age := age - init.year + as.integer(factor_to_char(year))] # align age with year
      output <- output[between(age, ageL, ageH)]
      output <- agegroup_aggregation_RF(output)
      
    }
    output[, population := NULL] # otherwise need to be divided by pop.fraction
    return(output)
  })

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

agegroup_aggregation <-
  cmpfun(function(x) {
    breaks <- c(ageL, 50L, 70L, ageH)
    labels <- c(paste0(ageL, "-49"), "50-69", paste0("70-", ageH))
    ff <- data.table(age = ageL:ageH)
    ff[, agegroup := cut(
      age,
      breaks = breaks,
      labels = labels,
      include.lowest = T,
      right = F,
      ordered_result = T)]
    
    x[ff, agegroup := i.agegroup, on = "age"]
    x[, age := NULL]
    # replace non finite
    for (j in outersect(names(x), c("agegroup", "year")))
      set(x,which(!is.finite(x[[j]])),j,0)
    
    
    x <- x[, lapply(.SD, sum, na.rm = TRUE), 
           keyby = .(year, agegroup)]
    return(x)
  }
  )

agegroup_aggregation_RF <-
  cmpfun(function(x) {
    breaks <- c(ageL, 50L, 70L, ageH)
    labels <- c(paste0(ageL, "-49"), "50-69", paste0("70-", ageH))
    ff <- data.table(age = ageL:ageH)
    ff[, agegroup := cut(
      age,
      breaks = breaks,
      labels = labels,
      include.lowest = T,
      right = F,
      ordered_result = T)]
    
    x[ff, agegroup := i.agegroup, on = "age"]
    x[, age := NULL]
    # replace non finite
    for (j in outersect(names(x), c("agegroup", "year", "scenario", "mc")))
      set(x,which(!is.finite(x[[j]])),j,0)
    
    
    x <- x[, lapply(.SD, function(x, y){ # weighted mean
      sum(x*y, na.rm = T) / sum(y, na.rm = T)
    }, population), 
           keyby = .(year, scenario, agegroup, mc)]
    return(x)  }
  )

