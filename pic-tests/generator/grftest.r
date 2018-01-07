library(RandomFields)
library(R.matlab)
library(psych)
rm(list=ls())
setwd("~/repo/convex-optimization/pic-tests/generator")

params = list(c(1, 0.25, 0.05), c(1, 1, 0.15), c(1, 2.5, 0.3), c(1, 0.5, 0.4), c(4, 4.5, 0.1))
resolutions = c(8, 16, 24, 32, 64, 96)
# resolutions = c(32)
index = 2
for (param in params) {
  picList = list()
  count = 1
  for (res in resolutions) {
    cat(sprintf('Res %d with param %f.\n', res, param[2]))
    mean = 0
    variance = param[0]
    nugget  = param[1]
    scale = 1/param[2]
    x = seq(0, 1, len=res)
    y = seq(0, 1, len=res)
    model = RMmatern(nu=param[1], scale=param[2])
    f = RFsimulate(model=model, x, y, n = 5)
    data1 = f$variable1.n1
    data2 = f$variable1.n2
    data3 = f$variable1.n3
    data4 = f$variable1.n4
    data5 = f$variable1.n5
    if (index == 5) {
      data1 = exp(data1)
      data2 = exp(data2)
      data3 = exp(data3)
      data4 = exp(data4)
      data5 = exp(data5)
    }
    if (index == 6) {
      data1 = logistic(data1*2)
      data2 = logistic(data2*2)
      data3 = logistic(data3*2)
      data4 = logistic(data4*2)
      data5 = logistic(data5*2)
    }
    data1 = data1 - min(data1) + 1e-10
    data1 = data1 / sum(data1)
    data2 = data2 - min(data2) + 1e-10
    data2 = data2 / sum(data2)
    data3 = data3 - min(data3) + 1e-10
    data3 = data3 / sum(data3)
    data4 = data4 - min(data4) + 1e-10
    data4 = data4 / sum(data4)
    data5 = data5 - min(data5) + 1e-10
    data5 = data5 / sum(data5)
    dim(data1) = c(res, res)
    dim(data2) = c(res, res)
    dim(data3) = c(res, res)
    dim(data4) = c(res, res)
    dim(data5) = c(res, res)
    picList[[count]] = as.matrix(rbind(data1, data2, data3, data4, data5))
    image(x, y, data1)
    count = count + 1
  }
  writeMat(paste(c("GRF_", index, ".mat"), collapse=""), picx8=picList[[1]], picx16=picList[[2]], picx24=picList[[3]], picx32=picList[[4]], picx64=picList[[5]], picx96=picList[[6]])
  index = index + 1
}

picList = list()
count = 1
for (res in resolutions) {
  cat(sprintf('Res %d.\n', res))
  x = seq(0, 1, len=res)
  y = seq(0, 1, len=res)
  model = RMcauchy(gamma=10)
  f = RFsimulate(model=model, x, y, n = 5)
  data1 = f$variable1.n1
  data2 = f$variable1.n2
  data3 = f$variable1.n3
  data4 = f$variable1.n4
  data5 = f$variable1.n5
  data1 = data1 - min(data1) + 1e-10
  data1 = data1 / sum(data1)
  data2 = data2 - min(data2) + 1e-10
  data2 = data2 / sum(data2)
  data3 = data3 - min(data3) + 1e-10
  data3 = data3 / sum(data3)
  data4 = data4 - min(data4) + 1e-10
  data4 = data4 / sum(data4)
  data5 = data5 - min(data5) + 1e-10
  data5 = data5 / sum(data5)
  dim(data1) = c(res, res)
  dim(data2) = c(res, res)
  dim(data3) = c(res, res)
  dim(data4) = c(res, res)
  dim(data5) = c(res, res)
  picList[[count]] = as.matrix(rbind(data1, data2, data3, data4, data5))
  image(x, y, data1)
  count = count + 1
}
writeMat(paste(c("GRF_", 7, ".mat"), collapse=""), picx8=picList[[1]], picx16=picList[[2]], picx24=picList[[3]], picx32=picList[[4]], picx64=picList[[5]], picx96=picList[[6]])
