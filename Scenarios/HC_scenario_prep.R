# predict undiagnosed diabetes
set(POP, NULL, "undiag.diab", 0L)
POP[diabtotr == "2", undiag.diab := pred.undiag.diab(.N, qimd)]

#QRisk score
set(POP, NULL, "qrisk2", 0)
set(POP, NULL, "diab_type_I", rbinom(nrow(POP), 1, 0.005))# diab type 1

POP[,
    qrisk2 := QRisk(
      age,
      sex,
      af,
      ra,
      kiddiag,
      bpmed,
      diab_type_I,
      diabtotr,
      bmival,
      origin,
      famcvd,
      tctohdl,
      omsysval,
      smoke_cat,
      townsend
    )]
POP[, qrisk_grp := cut(qrisk2, 
                       breaks = c(0, 0.1, 0.2, Inf), 
                       labels = c("qriskupto10", "qrisk10to20", "qrisk20+"), 
                       include.lowest = T, 
                       right = F)]
POP[between(age, 40, 74), agerange :=  cut(age, 
                                           breaks = c(40, 50, 60, 70, 75), 
                                           labels = c("40-49", "50-59", "60-69", "70-74"), 
                                           include.lowest = T, 
                                           right = F, 
                                           ordered_result = T)]
