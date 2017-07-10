#cmpfile("./risk factor trajectories.R")
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

load(file = "./Lagtimes/bmi.svylm.rda")
load(file = "./Lagtimes/chol.svylm.rda")
load(file = "./Lagtimes/sbp.svylm.rda")
load(file = "./Lagtimes/smok.active.svylr.rda")
load(file = "./Lagtimes/smok.exactive.svylr.rda")
load(file = "./Lagtimes/smok.cess.svylr.rda")
load(file = "./Lagtimes/smok.cess.success.log.rda")
load(file = "./Lagtimes/smok.cess.success.parabola.rda")
load(file = "./Lagtimes/smok.start.svylr.rda")
load(file = "./Lagtimes/cigdyal.svylr.rda")
load(file = "./Lagtimes/fv.svylr.rda")
load(file = "./Lagtimes/pa.svylr.rda")
load(file = "./Lagtimes/salt.rq.rda")
load(file = "./Lagtimes/tctohdl.svylm.rda")
load(file = "./Lagtimes/famcvd.svylr.rda")
load(file = "./Lagtimes/af.svylr.rda")
load(file = "./Lagtimes/kidfailgp.svylr.rda")
load(file = "./Lagtimes/bpmed.svylr.rda")
load(file = "./Lagtimes/undiag.diab.svylr.rda")
 
cat("Load RF trajectories\n")


# PA prediction ---------------------------------------
pred.pa <-
  cmpfun(
    function(year, age, sex, qimd, percentil.rank, fortune = 0.1) {
      percentil.rank <- jitter.constr(percentil.rank, fortune)
      
      newdata <- data.frame(
        year = year,
        age  = bound(age, 20, 84), 
        sex  = sex, 
        qimd = qimd
      )
      #code adapted from method getAnywhere(predict.polr)
      Terms <- delete.response(pa.svylr$terms)
      m <- model.frame(Terms, newdata, na.action = function(x) x, 
                       xlev = pa.svylr$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses"))) 
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts = pa.svylr$contrasts)
      xint <- match("(Intercept)", colnames(X), nomatch = 0L)
      if (xint > 0L)
        X <- X[, -xint, drop = FALSE]
      n <- nrow(X)
      q <- length(pa.svylr$zeta)
      eta <- drop(X %*% pa.svylr$coefficients)
      cc <- data.table(plogis(matrix(pa.svylr$zeta, n, q, byrow = TRUE) -
                                eta))
      set(cc, NULL, "V8", 1)
      #if (paired) set.seed(seed[[counter[[iterations]]]] + year)
      set(cc, NULL, "d", percentil.rank)
      for (k in 1:8) set(cc, NULL, k, cc[, Reduce(`<`, .SD), .SDcol = c(k, 9)])
      out <- cc[, Reduce(`+`, .SD), .SDcol = 1:8]

      # out <- predict_svylr(X, xint, pa.svylr$coefficients,
      #                      pa.svylr$zeta, percentil.rank )
      return(out) 
    } 
  )

# F&V prediction --------------------------------------
pred.fv <- 
  cmpfun(
    function(year, age, sex, qimd, percentil.rank, fortune = 0.1) {
      percentil.rank <- jitter.constr(percentil.rank, fortune) 
      
      newdata <- data.frame(
        year = year,
        age  = bound(age, 20, 84), 
        sex  = sex, 
        qimd = qimd
      )
      #code adapted from method getAnywhere(predict.polr)
      Terms <- delete.response(fv.svylr$terms)
      m <- model.frame(Terms, newdata, na.action = function(x) x, 
                       xlev = fv.svylr$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses"))) 
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts = fv.svylr$contrasts)
      xint <- match("(Intercept)", colnames(X), nomatch = 0L)
      if (xint > 0L)
        X <- X[, -xint, drop = FALSE]
      n <- nrow(X)
      q <- length(fv.svylr$zeta)
      eta <- drop(X %*% fv.svylr$coefficients)
      cc <- data.table(plogis(matrix(fv.svylr$zeta, n, q, byrow = TRUE) -
                                   eta))
      set(cc, NULL, "V9", 1)
      #if (paired) set.seed(seed[[counter[[iterations]]]] + year)
      set(cc, NULL, "d", percentil.rank)
      for (k in 1:9) set(cc, NULL, k, cc[, Reduce(`<`, .SD), .SDcol = c(k, 10)])
      out <- cc[, Reduce(`+`, .SD), .SDcol = 1:9]

      # crashes within foreach with error message "unimplemented type 'NULL' in 'coerceToInteger'"
      # out <- predict_svylr(X, xint, fv.svylr$coefficients, fv.svylr$zeta, percentil.rank)
     return(out) 
    }
  ) 


# test
# summary(factor(pred.fv(sample(c(0:50), n, replace = T), 
#                        sample(c(20,85), n, replace = T), 
#                        sample(c(1,2), n, replace = T), 
#                        sample(c(1:5), n, replace = T),
#                        runif(n, 10, 90),
#                        sample(c(1,10), n, replace = T))))/n
# SPOP2011[age>19, summary(factor(porftvg))/.N]


# Smoke initiation --------------------------------------------------------
# Gives the annual probability of a never smoker to become  smoker next year
# all other remain never smokers
pred.nev0sm1 <-
  cmpfun(
    function(year, age, sex) {
      
      pnev0sm1 <- 
        data.frame(
          predict(smok.start.svylr, 
                          data.table(year = year,
                                     age  = bound(age, 16, 84),
                                     sex  = sex),
                          type = "response", 
                          se.fit = T)
        )
      return(pnev0sm1[[1]])
    }
  )

# plot(pred.nev0sm1(0, 16:60, 1), ylim=c(0,0.5))
# lines(pred.nev0sm1(16:60, "3"), ylim=c(0,0.2))


# Smoke cessation ---------------------------------------------------------
# Predicts the annual probability of a smoker to become ex-smoker
pred.sm0ex1 <-
  cmpfun(
    function(year, age, sex, qimd) {
      qimd <- mapvalues(qimd, 1:5, c(1,2,2,2,3))
      sm0ex1 <- 
        data.frame(
          predict(smok.cess.svylr,
                        data.table(year = year,
                                   age = bound(age, 16, 84), 
                                   sex = sex, 
                                   qimd = qimd), 
                        type = "response",
                        se.fit = T)
        )
      return(sm0ex1[[1]])
    }
  )

# for (jj in 1:5) {
#     plot(16:90, pred.sm0ex1(0, 16:90, 1, jj), ylim=c(0,0.4))
# }


# Smoke relapse -----------------------------------------------------------
# # Predicts probability of ex-smoker to become active smoker (relapse) (only works for 1<endsmoke<10). Else should be 0
pred.ex0sm1 <- cmpfun(
  function(endsmoke, sex, qimd, type = c("parabola", "log")) {
    if (type == "parabola") {
      ex0sm1 <-
        data.table(
          predict(
            smok.cess.success.parabola,
            data.table(endsmoke = bound(endsmoke, 1, 20),
                       sex = sex,
                       qimd = qimd),
            type = "response",
            se.fit = F))
    }
    if (type == "log") {
      ex0sm1 <-
        data.table(
          predict(
            smok.cess.success.log,
            data.table(endsmoke = bound(endsmoke, 1, 20),
                       sex = sex,
                       qimd = qimd),
            type = "response",
            se.fit = F))
    }
    return(ex0sm1)
  }
)
#pred.ex0sm1(1:10, 1, 1, "parabola")

# Smoke prevalence --------------------------------------------------------
# predicts the active smoker prevalence
# fitted for ages 16 to 20 inclusive
pred.sm0prev <- cmpfun(function(year, age, sex, qimd) {
  sm0prev <- predict(smok.active.svylr,
                     data.table(year = year, 
                                age = bound(age, 16, 84), 
                                sex = sex, 
                                qimd = qimd),
                     type = "response",
                     se.fit=F)
  #return(sm0prev[[1]])
  return(rbinom(length(sm0prev), 1, sm0prev))
}
)

# Ex-Smoke prevalence --------------------------------------------------------
# predicts the ex smoker prevalence
# fitted for ages 16 to 20 inclusive
pred.exsm0prev <- cmpfun(function(year, age, sex, qimd) {
  sm0prev <- predict(smok.exactive.svylr,
                     data.table(year = year, 
                                age = bound(age, 16, 84), 
                                sex = sex, 
                                qimd = qimd),
                     type = "response",
                     se.fit=F)
  #return(sm0prev[[1]])
  return(rbinom(length(sm0prev), 1, sm0prev))
}
)

# cigdyal prediction ---------------------------------------
pred.cigdyal <-
  cmpfun(
    function(year, age, sex, qimd, smokyrs, percentil.rank, fortune = 0.1) {
      percentil.rank <- jitter.constr(percentil.rank, fortune)
      
      newdata <- data.frame(
        year = year,
        age  = bound(age, 16, 84), 
        sex  = sex, 
        qimd = qimd,
        smokyrs = bound(smokyrs, 0, 76)
      )
      #code adapted from method getAnywhere(predict.polr)
      Terms <- delete.response(cigdyal.svylr$terms)
      m <- model.frame(Terms, newdata, na.action = function(x) x, 
                       xlev = cigdyal.svylr$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses"))) 
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts = cigdyal.svylr$contrasts)
      xint <- match("(Intercept)", colnames(X), nomatch = 0L)
      # if (xint > 0L) 
      #   X <- X[, -xint, drop = FALSE]
      # n <- nrow(X)
      # q <- length(cigdyal.svylr$zeta)
      # eta <- drop(X %*% cigdyal.svylr$coefficients)
      # cc <- data.table(plogis(matrix(cigdyal.svylr$zeta, n, q, byrow = TRUE) - 
      #                              eta))
      # set(cc, NULL, "V4", 1)
      # set(cc, NULL, "d", percentil.rank)
      # for (k in 1:4) set(cc, NULL, k, cc[, Reduce(`<`, .SD), .SDcol = c(k, 5)])
      # cc[, out := Reduce(`+`, .SD), .SDcol = 1:4]
      # cc[, out := out * 10]
      # cc[out == 0, out := 5]
      # return(cc$out) 
      
      return(
        bound(
        10 * predict_svylr(X, xint, cigdyal.svylr$coefficients,
                           cigdyal.svylr$zeta, percentil.rank),
        5, Inf, T)
      )
    } 
  )


# BMI prediction ----------------------------------------------------------
# Define function for bmi projection (predicts mean bmi)
pred.bmi <- 
  cmpfun(
    function(year, age, sex, qimd, a30to06m) {
      
      pr <- data.frame(
        predict(
          bmi.svylm, 
          data.frame(
            year         = year,
            age          = bound(age, 20, 84), 
            sex          = sex, 
            qimd         = qimd,
            a30to06m.imp = a30to06m
          ), 
          type = "response", 
          se.fit = T
        )
      )
      return(pr[[1]])
    }
  )

# test
# summary(pred.bmi(sample(c(0:50), n, replace = T), 
#                  sample(c(20,85), n, replace = T), 
#                  sample(c(1,2), n, replace = T), 
#                  sample(c(1:5), n, replace = T),
#                  sample(c(1,10), n, replace = T)))


# Salt prediction ---------------------------------------------------------
# Returns a dataframe of 24h salt percentiles by year, age, sex, qimd 
pred.salt <- 
  cmpfun(
    function(year, lag = cancer.lag) {
      year <- 
        switch(EXPR = as.character((year - lag)),
               "-16" =  lag - 8,
               "-15" =  lag - 8.5, 
               "-14" =  lag - 9,
               "-13" =  lag - 9.5, 
               "-12" =  lag - 9, #9
               "-11" =  lag - 8.5, #8.5
               "-10" =  lag - 8, #8
               "-9"  =  lag - 7.7,
               "-8"  =  lag - 7.5,
               year
        )
      
      tmp <- expand.grid(
        year = year-lag,
        age  = (19 - lag):(ageH - lag),
        #age  = (ageL-lag):(ageH-lag),
        sex  = factor(1:2),
        qimd = ordered(1:5)
      )
      cc <- predict(salt.rq, tmp)^3
      
      tmp <- data.table(cbind(tmp, cc))
      tmp[, `:=` (year = NULL, age = age + lag)]
      setnames(tmp,
               paste0("tau= ", sprintf("%.2f", c(0.01, 1:19/20, 0.99))),
               paste0(c(0.01, 1:19/20, 0.99)))
      
      tmp <- melt(tmp, 1:3, 
                  variable.name = "percentile",
                  value.name = "salt.u", 
                  variable.factor = F)
      
      tmp[, percentile := as.numeric(percentile)]
      tmp[, salt.l := shift(salt.u, 1, 1, "lag"),
          by = .(age, sex, qimd)]
      tmp[salt.u<salt.l, salt.t := salt.u] # logic to reverse column l u columns
      tmp[salt.u<salt.l, `:=` (salt.u = salt.l, salt.l = salt.t)]
      if ("salt.t" %in% names(tmp)) tmp[, salt.t := NULL]
      return(tmp)
    }
  )


# SBP prediction ----------------------------------------------------------
# Define function for sbp projection (for DT needs the by= to work correctly with mean(bmival)) (predicts mean sbp)
pred.sbp <-
  cmpfun(
    function(year, age, sex, qimd, bmival, cigst1, a30to06m) {
      # cigst2 <- factor(as.integer(cigst1 == "4"))
      cigst2 <- mapvalues(cigst1,  c(4:1 ), c(1,0,0,0))
      bmival <- bound(bmival, 16, 40) # otherwise predicts NAN values
      pr <- data.frame(
        predict(
          sbp.svylm, 
          data.table(
            year = year, 
            age = bound(age, 20, 84), 
            sex = sex,
            qimd = qimd, 
            bmival = bmival,
            cigst2 = cigst2,
            a30to06m.imp = a30to06m),
          type = "response", se.fit = T
        )
      )
      return(pr[[1]])
    }
  )

#test
# summary(pred.sbp(sample(c(0:50), n, replace = T), 
#                  sample(c(20,85), n, replace = T), 
#                  sample(c(1,2), n, replace = T), 
#                  sample(c(1:5), n, replace = T),
#                  runif(n, 10, 90),
#                  sample(c(1,10), n, replace = T)))


# Chol prediction ---------------------------------------------------------
# Define function for chol projection (for ages above 30)
pred.chol <- cmpfun(function(year, age, sex, qimd, bmival) {
  bmival <- bound(bmival, 19, 35) # otherwise predicts NAN values
  pr <- data.frame(
    predict(
      chol.svylm, 
      data.table(
        year = year,
        age = bound(age, 20, 84), 
        sex = sex, 
        qimd = qimd, 
        bmival = bmival),
      type = "response", 
      se.fit=T
    )
  )
  return(pr[[1]])
  #return(rtruncnorm(nrow(pr), a = 2.5, b = 12,  pr[[1]], pr[[2]]))
  #return(rnorm(nrow(pr), pr[[1]], pr[[2]]))
}
)

# test
# summary(pred.chol(sample(c(-10:50), n, replace = T), 
#                   sample(c(20,85), n, replace = T), 
#                   sample(c(1,2), n, replace = T), 
#                   sample(c(1:5), n, replace = T),
#                   runif(n, 10, 50),runif(n, 10, 20),
#                   sample(c(1,10), n, replace = T)))

# TC to HDL prediction ---------------------------------------------------------
# Define function for hdl estimation
pred.tctohdl <- 
  cmpfun(
    function(cholval1, age, sex, qimd, bmival, a30to06m, cigst1) {
  bmival <- bound(bmival, 16, 40) # otherwise predicts NAN values
  cigst2 <- mapvalues(cigst1,  c(4:1 ), c(1,0,0,0))
  pr <- data.frame(
    predict(
      tctohdl.svylm, 
      data.frame(
        cholval1     = bound(cholval1, 2, 12),
        age          = bound(age, 20, 84), 
        sex          = sex, 
        qimd         = qimd, 
        bmival       = bmival,
        a30to06m.imp = a30to06m,
        cigst2       = cigst2),
      type = "response", 
      se.fit = T
    )
  )
  return(pr[[1]])
  #return(rtruncnorm(nrow(pr), a = 2.5, b = 12,  pr[[1]], pr[[2]]))
  #return(rnorm(nrow(pr), pr[[1]], pr[[2]]))
}
)

# FamCVD prediction ---------------------------------------------
pred.famcvd <- cmpfun(function(n, age, qimd) {
  newdata <-
    data.table(
      age          = bound(age, 20, 84), 
      qimd         = qimd
    )
  
  type <-"response"
  total <- NULL
  tt <- delete.response(terms(formula(famcvd.svylr)))
  mf <- model.frame(tt, data = newdata)
  mm <- model.matrix(tt, mf)
  if (!is.null(total) && attr(tt, "intercept")) {
    mm[, attr(tt, "intercept")] <- mm[, attr(tt, "intercept")] * 
      total
  }
  eta <- drop(mm %*% coef(famcvd.svylr))
  eta <- switch(type, link = eta, response = famcvd.svylr$family$linkinv(eta))
  rbinom(n, 1, eta)
  #return(rtruncnorm(nrow(pr), a = 0, b = 1, mean=pr[[1]], sd=pr[[2]])) 
}
)

# AF prevalence prediction ---------------------------------------------
pred.af <- cmpfun(function(n, age, omsysval) {
  newdata <-
    data.table(
      age          = bound(age, 20, 84), 
      omsysval     = omsysval
    )
  
  type <-"response"
  total <- NULL
  tt <- delete.response(terms(formula(af.svylr)))
  mf <- model.frame(tt, data = newdata)
  mm <- model.matrix(tt, mf)
  if (!is.null(total) && attr(tt, "intercept")) {
    mm[, attr(tt, "intercept")] <- mm[, attr(tt, "intercept")] * 
      total
  }
  eta <- drop(mm %*% coef(af.svylr))
  eta <- switch(type, link = eta, response = af.svylr$family$linkinv(eta))
  rbinom(n, 1, eta)
  #return(rtruncnorm(nrow(pr), a = 0, b = 1, mean=pr[[1]], sd=pr[[2]])) 
}
)

# Kidney disease prevalence prediction ---------------------------------------------
pred.kidfailgp <- cmpfun(function(n, age, omsysval) {
  newdata <-
    data.table(
      age          = bound(age, 20, 84), 
      omsysval     = omsysval    )
  
  type <-"response"
  total <- NULL
  tt <- delete.response(terms(formula(kidfailgp.svylr)))
  mf <- model.frame(tt, data = newdata)
  mm <- model.matrix(tt, mf)
  if (!is.null(total) && attr(tt, "intercept")) {
    mm[, attr(tt, "intercept")] <- mm[, attr(tt, "intercept")] * 
      total
  }
  eta <- drop(mm %*% coef(kidfailgp.svylr))
  eta <- switch(type, link = eta, response = kidfailgp.svylr$family$linkinv(eta))
  rbinom(n, 1, eta)
  #return(rtruncnorm(nrow(pr), a = 0, b = 1, mean=pr[[1]], sd=pr[[2]])) 
}
)

# BP medication prediction ---------------------------------------------
pred.bpmed <- cmpfun(function(n, year, age, sex, qimd, omsysval) {
  newdata <-
    data.frame(
      year         = year,
      age          = bound(age, 20, 84), 
      sex          = sex,
      qimd         = qimd,
      omsysval     = bound(omsysval, 70, 220)
    )
  
  type <-"response"
  total <- NULL
  tt <- delete.response(terms(formula(bpmed.svylr)))
  mf <- model.frame(tt, data = newdata)
  mm <- model.matrix(tt, mf)
  if (!is.null(total) && attr(tt, "intercept")) {
    mm[, attr(tt, "intercept")] <- mm[, attr(tt, "intercept")] * 
      total
  }
  eta <- drop(mm %*% coef(bpmed.svylr))
  eta <- switch(type, link = eta, response = bpmed.svylr$family$linkinv(eta))
  rbinom(n, 1, eta)
  #return(rtruncnorm(nrow(pr), a = 0, b = 1, mean=pr[[1]], sd=pr[[2]])) 
}
)

# Undiagnosed  prediction ---------------------------------------------
pred.undiag.diab <- 
  cmpfun(
    function(n, qimd) {
  newdata <-
    data.frame(
      qimd = qimd)
  
  type <-"response"
  total <- NULL
  tt <- delete.response(terms(formula(undiag.diab.svylr)))
  mf <- model.frame(tt, data = newdata)
  mm <- model.matrix(tt, mf)
  if (!is.null(total) && attr(tt, "intercept")) {
    mm[, attr(tt, "intercept")] <- mm[, attr(tt, "intercept")] * 
      total
  }
  eta <- drop(mm %*% coef(undiag.diab.svylr))
  eta <- switch(type, link = eta, response = undiag.diab.svylr$family$linkinv(eta))
  rbinom(n, 1, eta)
  #return(rtruncnorm(nrow(pr), a = 0, b = 1, mean=pr[[1]], sd=pr[[2]])) 
}
)

