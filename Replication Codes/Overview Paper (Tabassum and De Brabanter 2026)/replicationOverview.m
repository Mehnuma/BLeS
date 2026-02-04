
data = readmatrix('simEx1.csv');
qdata = readmatrix('qregData.csv');
pdata = readmatrix('pData.csv');
nsims = 100;
n = 200;

%% HHJ Method
time_hhj = nan(nsims,1);
index_hhj = 1;
ln_hhj = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_hhj:(index_hhj+n-1),1);
    ln_hhj(iter) = HHJ_method(y,3,30,'mean','variance',500,0,0,[]); %l_init =5, m = 30, B = 800; %BLeS
    index_hhj = index_hhj+n;
    time_hhj(iter) = toc;
end

%% BK Method
time_bk = nan(nsims,1);
index_bk = 1;
ln_bk = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_bk:(index_bk+n-1),1);
    ln_bk(iter) = BK_method(y,'median'); %BLeS
    index_bk = index_bk+n;
    time_bk(iter) = toc;
end

%% cPW Method (Circular BB)
time_cbb_patton = nan(nsims,1);
index_cbb_patton = 1;
ln_cbb_patton = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_cbb_patton:(index_cbb_patton+n-1),1);
    ln_cbb_patton(iter) = round(opt_BL_Patton(y)); %Patton
    index_cbb_patton = index_cbb_patton+n;
    time_cbb_patton(iter) = toc;
end

time_cbb_bles = nan(iter,1);
index_cbb = 1;
ln_cbb = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_cbb:(index_cbb+n-1),1);
    ln_cbb(iter) = cPW_method(y,'circularBB'); %BLeS
    index_cbb = index_cbb+n;
    time_cbb_bles(iter) = toc; 
end

%% cPW Method (Stationary BB)
time_sbb_patton = nan(nsims,1);
index_sbb_patton = 1;
ln_sbb_patton = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_sbb_patton:(index_sbb_patton+n-1),1);
    ln_sbb_patton(iter) = opt_BL_Patton(y); %Patton
    index_sbb_patton = index_sbb_patton+n;
    time_sbb_patton(iter) = toc;
end

time_sbb_bles = nan(nsims,1);
index_sbb = 1;
ln_sbb = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_sbb:(index_sbb+n-1),1);
    ln_sbb(iter) = cPW_method(y,'stationaryBB'); %BLeS
    index_sbb = index_sbb+n;
    time_sbb_bles(iter)=toc;
end

%% NPPI Method
time_nppi = nan(nsims,1);
index_nppi = 1;
ln_nppi = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_nppi:(index_nppi+n-1),1);
    ln_nppi(iter) = NPPI_method(y,'mean','quantile',500,0.8238,threshold,0.8); %BLeS
    index_nppi = index_nppi+n;
    time_nppi(iter) = toc;
end

%% Block length selection for TBB
time_tbb = nan(nsims,1);
index_tbb = 1;
ln_tbb = nan(nsims,1);
for iter=1:nsims
    tic
    y = data(index_tbb:(index_tbb+n-1),1);
    ln_tbb(iter) = TBB_bles(y,'TBB'); %BLeS
    index_tbb = index_tbb+n;
    time_tbb(iter) = toc;
end

%% Block length selection for SETBB
time_setbb = nan(nsims,1);
index_setbb = 1;
ln_setbb = nan(nsims,1);
for iter=1:nsims
    tic
    X = qdata(index_setbb:(index_setbb+n-1),1:2);    
    Y = qdata(index_setbb:(index_setbb+n-1),3);
    ln_setbb(iter) = SETBB_bles(X,Y,0.8,'SETBB',500); %BLeS
    index_setbb = index_setbb+n;
    time_setbb(iter) = toc;
end

%% Block length selection for the BD Method
time_bd = nan(nsims,1);
index_bd = 1;
ln_bd = nan(nsims,1);
for iter=1:nsims
    tic
    y = pdata(index_bd:(index_bd+n-1),1);
    ln_bd(iter) = BD_bles(y,6,'overall mean','EMBB'); %BLeS
    index_bd = index_bd+n;
    time_bd(iter) = toc; 
end
