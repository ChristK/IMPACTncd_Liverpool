load("./Scenarios/Tables.rda")
salt.reform.effect <- copy(Tables$sbp.SAQ[scenario %in% c("current trends", "stricter policy") & agegroup %in% c( "30-34", "35-39", 
                                                                                                                 "40-44", "45-49", "50-54",
                                                                                                                 "55-59", "60-64", "65-69",
                                                                                                                 "70-74", "75-79", "80-84", "85+")])
salt.reform.effect[, year := year.cvdlag]
salt.reform.effect[, year.cvdlag := NULL]
salt.reform.effect <- salt.reform.effect[between(year, 2016, 2035)]
# salt.reform.effect[, `:=` (
#   lui = lui / mean,
#   uui = uui / mean
# )]
tt <- split(salt.reform.effect, f = salt.reform.effect$scenario)
salt.reform.effect <- tt$`current trends`[tt$`stricter policy`, on = c("year", "agegroup", "sex", "qimd"), .(
  year, agegroup, sex, qimd, 
  mean = i.mean - mean,
  uui = i.uui - uui,
  lui = i.lui - lui
)]
salt.reform.effect[mean>0, mean := 0]
salt.reform.effect[lui>0, lui := 0]
salt.reform.effect[uui>0, uui := 0]
salt.reform.effect[, max := base::pmax(mean, lui, uui) ]
salt.reform.effect[, min := base::pmin(mean, lui, uui) ]
salt.reform.effect[, mid := mean + lui + uui - max - min]
salt.reform.effect[, c("mean", "uui", "lui") := NULL]


tt <- salt.reform.effect[agegroup %in% c("30-34")]
tt[, agegroup := "25-29"]
salt.reform.effect <- rbind(salt.reform.effect, tt)
tt[, agegroup := "20-24"]
salt.reform.effect <- rbind(salt.reform.effect, tt)
tt[, agegroup := "20-24"]
tt <- salt.reform.effect[year == 2035]
l <- rep(list(tt), 10)
tt <- rbindlist(l, idcol = "id")
tt[, year:= year + id - 1L]
tt[, id := NULL]
salt.reform.effect <- rbind(salt.reform.effect, tt)

salt.reform.effect[, agegroup := factor(agegroup)]
salt.reform.effect[sex == "Men", sex := "1"]
salt.reform.effect[sex == "Women", sex := "2"]

setkey(salt.reform.effect, year, agegroup, sex, qimd)
fwrite(salt.reform.effect, "./Scenarios/salt_reform_effect.csv")


