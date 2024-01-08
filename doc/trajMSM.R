## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(trajMSM)

## -----------------------------------------------------------------------------
Obsdata = gendatatraj()
Rdat = longformat(Obsdata = Obsdata, varying = 1:5)
head(Rdat)

## ----message=TRUE, warning=FALSE----------------------------------------------

res.traj = buildtraj(Rdat = Rdat, k=3,formula = cbind(A,1-A) ~ time, id="id")
dpost = res.traj$dpost
head(dpost)
plotraj(Rdat = Rdat, dpost = dpost, formula = A ~ time + class,
trt = "A",time = "time",id="id",class = "class", FUN = "mean")


