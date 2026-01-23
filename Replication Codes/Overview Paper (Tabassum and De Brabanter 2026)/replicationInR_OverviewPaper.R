
data = read.csv('simEx1.csv', header = F)
qData = read.csv('qregData.csv')
pData = read.csv('pData.csv')

nsim = 100
n = 200

library(blocklength)
library(np)
library(OBL)
library(QregBB)
library(boodd)

#****** Package: blocklength *****
index_hhj = 1
ln_hhj_blocklength = NULL
ln_cbb_blocklength = NULL
ln_sbb_blocklength = NULL
for(iter in (1:nsim)){
  y = data[index_hhj:(index_hhj+n-1),1]
  hhj1_out = blocklength::hhj(y, pilot_block_length =5, sub_sample = 30, k = "bias/variance", nb=500)
  ln_hhj_blocklength[iter] = hhj1_out$`Optimal Block Length`
  hhj2_out = blocklength::pwsd(y,correlogram=FALSE)
  ln_sbb_blocklength[iter] = hhj2_out$BlockLength[1]
  ln_cbb_blocklength[iter] = round(hhj2_out$BlockLength[2])
  index_hhj = index_hhj+n;
}

#****** Package: np *****
index_np = 1
ln_cbb_np = NULL
ln_sbb_np = NULL
for(iter in (1:nsim)){
  y = data[index_np:(index_np+n-1),1]
  np_out = np::b.star(y)
  ln_sbb_np[iter] = np_out[1]
  ln_cbb_np[iter] = round(np_out[2])
  index_np = index_np+n;
}

#****** Package: OBL *****
index_obl = 1
ln_cbb_obl = NULL
ln_tbb_obl = NULL
for(iter in (1:nsim)){
  y = data[index_obl:(index_obl+n-1),1]
  OBL_out = OBL::blockboot(y,R=2,seed=6,n_cores=1)
  ln_cbb_obl[iter] = OBL_out$lb[3]
  ln_tbb_obl[iter] = OBL_out$lb[4]
  index_obl = index_obl+n;
}

#****** Package: QregBB *****
index_qregbb = 1
ln_qregbb = NULL
for(iter in (1:nsim)){
  X = qData[index_qregbb:(index_qregbb+n-1),1:2]
  X = as.matrix(X)
  Y = qData[index_qregbb:(index_qregbb+n-1),3]
  qregbb_out = QregBB::getNPPIblksizesQR(Y,X,0.8)
  ln_qregbb[iter] = qregbb_out$l.opt.SETBB
  index_qregbb = index_qregbb+n;
}

#****** Package: boodd *****
index_boodd = 1
bopt_boodd = NULL
period = 12

for (iter in (1:nsim)){
  y = pData[index_boodd:(index_boodd+n-1),1]
  bopt_boodd[iter] =  boodd::bopt_circy(y, period, parameter= "mean", method= "GSBB")
  index_boodd = index_boodd+n
}
  

