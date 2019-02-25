# Borges_etal_analysis_congo_CCCexcluded.R 

# this script reads in width and length data from shapefiles
# and extrapolates them to make a new estimate of RSSA by 
# stream order in the congo basin (with the CCC removed):


library(foreign)

################################################################################
# specify shapefile dbf paths:
################################################################################
wd = "E:/research/2018_12_20_Borges_etal_Congo/git/"

hS_path = paste0(wd, "in/hSHEDS/hSHEDs_grwl_ccc_glwd.dbf")
hS_oDir = paste0(wd, "out/hSHEDS")

# if needed copy input shapefile files to output directory: 
source_dir = dirname(hS_path)
source_paths = list.files(source_dir, full.names=T)
targ_paths = sub(source_dir, hS_oDir, source_paths)
if (F %in% file.exists(targ_paths)){
  file.copy(source_paths, targ_paths)
}
hS_oPath = targ_paths[grep(".dbf", targ_paths)]

# read in shapefile attribute table: 
hS = read.dbf(hS_path) # hSHEDS flowlines with GRWL & GLCD2009 data joined


# Calculate area of congo basin minus area of CCC from 
# shapefile (area calculated using UTM 34N projection):
ccc = read.dbf("E:/research/2018_12_20_Borges_etal_Congo/GIS/output/ccc/glc2009_ccc_v1.dbf")
outsideoutsideCCC_km2 = 3689190-sum(ccc$area)


# increase stream order by one to make it consistent with Raymond et al:
hS$ORDER_ = hS$ORDER_ + 1

################################################################################
# river width by stream order:
################################################################################


# hSHEDS with stream order > 3, GRWL widths > 0, non lakes, & outside the CCC:
hS_OrdGT3 = hS[hS$ORDER_>4 & 
                 hS$width_med_>0 &
                 hS$GLWD_ID == 0 &
                 hS$gridcode == 0, ] 


# widths derived from DHG:
boxplot(hS_OrdGT3$WIDTH ~ hS_OrdGT3$ORDER_, 
        ylim = range(hS_OrdGT3$width_med_),
        main="Outside CCC SRTM-DHG",
        xlab="Order",
        ylab="Width (m)",
        col=rgb(1,0,0,0.3))
# widths measured by GRWL:
boxplot(hS_OrdGT3$width_med_ ~ hS_OrdGT3$ORDER_, 
        main="Outside CCC GRWL",
        xlab="Order",
        ylab="Width (m)",
        col=rgb(0,0,0,0.3))

# median GRWL width by order: 
uniqOrd = unique(hS_OrdGT3$ORDER_)
uniqOrd = sort(uniqOrd)

width_byOrd_obs = rep(NA, length(uniqOrd))
for (i in 1:length(uniqOrd)){
  width_byOrd_obs[i] = median(hS_OrdGT3$width_med_[hS_OrdGT3$ORDER_==uniqOrd[i]])
}

barplot(width_byOrd_obs, names.arg=uniqOrd, 
        main="Outside CCC Median Observed River Width by Order",
        xlab="Order",
        ylab="Median River Width (m)")

# extrapolate median widths down to 1st order streams using an expoential scaling:
logW_lm = lm(log(width_byOrd_obs) ~ uniqOrd)

# determine which orders need to be modeled:
modOrds = c(1:c(min(uniqOrd)-1))
modW = exp(logW_lm$coefficients[[2]]*modOrds+logW_lm$coefficients[[1]])
width_byOrd = c(modW, width_byOrd_obs)

# median widths according to the statistical stream order model:
modOrds_all = 1:max(uniqOrd)
modWidth_byOrd = exp(logW_lm$coefficients[[2]]*modOrds_all+logW_lm$coefficients[[1]])

# plot width extrapolation:
plot(range(c(modOrds, uniqOrd)), range(c(width_byOrd, modWidth_byOrd)), type='n',
     main="Outside CCC Median Width by Order",
     xlab="Order",
     ylab="Median Width (m)")
xSeq = seq(-1, max(uniqOrd)+1, length.out=100)
ySeq = exp(logW_lm$coefficients[[2]]*xSeq+logW_lm$coefficients[[1]])
lines(xSeq, ySeq, col=2, lty=1)
points(modOrds, modW, pch=16, col=2)
points(uniqOrd, width_byOrd_obs, pch=16, col=1)
legend("topleft", c("Observed", "Modeled"), col=c(1,2), pch=16, inset=0.05)


################################################################################
# river length by stream order:
################################################################################
# calculate median length (km) by order from hSHEDS flowlines:
uniqOrd = unique(hS$ORDER_)
uniqOrd = sort(uniqOrd)


# remove lakes and non CCC regions from lengths:
hS_noLakes = hS[hS$GLWD_ID == 0 & 
                  hS$gridcode == 0, ]

len_byOrd_obs = rep(NA, length(uniqOrd))
for (i in 1:length(uniqOrd)){
  len_byOrd_obs[i] = sum(hS_noLakes$LENGTH_KM[hS_noLakes$ORDER_==uniqOrd[i]] )
}

barplot(len_byOrd_obs, names.arg=uniqOrd, 
        main="Outside CCC Sum River Length",
        xlab="Order",
        ylab="Sum River Length (km)")

# extrapolate length down to 0-order streams: 
logL_lm = lm(log(len_byOrd_obs) ~ uniqOrd)

# determine which orders need to be modeled:
modOrds = c(1:c(min(uniqOrd)-1))
modL = exp(logL_lm$coefficients[[2]]*modOrds+logL_lm$coefficients[[1]])
len_byOrd = c(modL, len_byOrd_obs)

# median length according to the statistical stream order model:
modOrds_all = 1:max(uniqOrd)
modLen_byOrd = exp(logL_lm$coefficients[[2]]*modOrds_all+logL_lm$coefficients[[1]])

# plot width extrapolation:
plot(c(modOrds, uniqOrd), len_byOrd, type='n',
     main="Outside CCC Median Length by Order",
     xlab="Order",
     ylab="Median Length (km)")
xSeq = seq(-1, max(uniqOrd)+1, length.out=100)
ySeq = exp(logL_lm$coefficients[[2]]*xSeq+logL_lm$coefficients[[1]])
lines(xSeq, ySeq, col=4)
points(modOrds, modL, pch=16, col=4)
points(uniqOrd, len_byOrd_obs, pch=16, col=1)
legend("topright", c("Observed", "Modeled"), col=c(1,4), pch=16, inset=0.05)



################################################################################
# river surface area by stream order:
################################################################################
RSSA_byOrd = len_byOrd * (width_byOrd*1e-3) # convert width from m to km


barplot(RSSA_byOrd, names.arg=modOrds_all, 
        main="Outside CCC River & Stream Surface Area",
        xlab="Order",
        ylab="River & Stream Surface Area (sq km)")

# create a RSSA by order table:
RSSA_byOrd_tab = as.data.frame(cbind(modOrds_all, round(width_byOrd, 1), round(len_byOrd), round(RSSA_byOrd)))
names(RSSA_byOrd_tab) = c("Order", "Med_Width_m", "Sum_Len_km", "Sum_Area_km2")
write.csv(RSSA_byOrd_tab, "E:/research/2018_12_20_Borges_etal_Congo/tabs/RSSA_byOrder_congo_CCCexcluded.csv", row.names=F)

# print basin %RSSA:
print(paste0(round(100*sum(RSSA_byOrd)/outsideCCC_km2, 2), "% RSSA in the outside the CCC"))




################################################################################
# calculate median slope by order:
################################################################################
# calculate median slope (grad) by order from hSHEDS flowlines:
uniqOrd = unique(hS$ORDER_)
uniqOrd = sort(uniqOrd)

# remove lakes and non CCC regions from lengths:
hS_noLakes = hS[hS$GLWD_ID == 0 & 
                  hS$gridcode == 0, ]


slope_byOrd_obs = rep(NA, length(uniqOrd))
for (i in 1:length(uniqOrd)){
  slope_byOrd_obs[i] = median(hS_noLakes$SLOPE[hS_noLakes$ORDER_==uniqOrd[i]] )
}

barplot(slope_byOrd_obs, names.arg=uniqOrd, 
        main="Outside CCC Median River Slope",
        xlab="Order",
        ylab="Median River Slope (grad)")

# extrapolate down to 0-order streams using a POWER LAW:
logS_lm = lm(log(slope_byOrd_obs) ~ log(uniqOrd))

# determine which orders need to be modeled:
modOrds = c(1:c(min(uniqOrd)-1))
modS = exp(logS_lm$coefficients[[1]])*modOrds^logS_lm$coefficients[[2]]
slope_byOrd = c(modS, slope_byOrd_obs)

# median  according to the statistical stream order model:
modOrds_all = 1:max(uniqOrd)
modSlope_byOrd = exp(logS_lm$coefficients[[1]])*modOrds_all^logS_lm$coefficients[[2]]

# plot extrapolation:
plot(c(modOrds, uniqOrd), slope_byOrd, type='n',
     main="Outside CCC Median River Slope by Order",
     xlab="Order",
     ylab="Median Slope (grad)", log="")
xSeq = seq(-1, max(uniqOrd)+1, length.out=100)
ySeq = exp(logS_lm$coefficients[[1]])*xSeq^logS_lm$coefficients[[2]]
lines(xSeq, ySeq, col=4)
points(modOrds, modS, pch=16, col=4)
points(uniqOrd, slope_byOrd_obs, pch=16, col=1)
legend("topright", c("Observed", "Modeled"), col=c(1,4), pch=16, inset=0.05)



################################################################################
# calculate median velocity by order:
################################################################################
# for hSHEDS with stream order > 3 and GRWL widths > 0 & non lakes, 
# use GRWL widths, otherwise use Andreadis et al's widths: 
OrdGT3 = hS$ORDER_>4 & hS$width_med_>0 & hS$GLWD_ID == 0

# calculate velocity (m/s) using Manning's formula:
hS$WIDTH[OrdGT3] = hS$width_med_[OrdGT3]

# hydraulic radius (m):
R = hS$WIDTH * hS$DEPTH / (2*hS$DEPTH + hS$WIDTH)

# manning's equation for flow velocity (m/s):
N = 0.035
hS$VELOCITY = N^(-1) * R^(2/3) * hS$SLOPE^(1/2)

# calculate median by order from hSHEDS flowlines:
uniqOrd = unique(hS$ORDER_)
uniqOrd = sort(uniqOrd)

# remove lakes and non CCC regions from lengths:
hS_noLakes = hS[hS$GLWD_ID == 0 & 
                  hS$gridcode == 0, ]

vel_byOrd_obs = rep(NA, length(uniqOrd))
for (i in 1:length(uniqOrd)){
  vel_byOrd_obs[i] = median(hS_noLakes$VELOCITY[hS_noLakes$ORDER_==uniqOrd[i]] )
}

barplot(vel_byOrd_obs, names.arg=uniqOrd, 
        main="Outside CCC Median Velocity",
        xlab="Order",
        ylab="Median River Velocity (m/s)")

# extrapolate down to 0-order streams using a exponetial fit:
logV_lm = lm(log(vel_byOrd_obs) ~ uniqOrd)

# determine which orders need to be modeled:
modOrds = c(1:c(min(uniqOrd)-1))
modV = exp(logV_lm$coefficients[[2]]*modOrds+logV_lm$coefficients[[1]])
vel_byOrd = c(modV, vel_byOrd_obs)

# median according to the statistical stream order model:
modOrds_all = 1:max(uniqOrd)
modVel_byOrd = exp(logV_lm$coefficients[[2]]*modOrds+logV_lm$coefficients[[1]])

# plot extrapolation:
plot(c(modOrdV, uniqOrd), vel_byOrd, type='n',
     main="Outside CCC Median River Velocity by Order",
     xlab="Order",
     ylab="Median Velocity (m/s)", log="")
xSeq = seq(-1, max(uniqOrd)+1, length.out=100)
ySeq = exp(logV_lm$coefficients[[2]]*xSeq+logV_lm$coefficients[[1]])
lines(xSeq, ySeq, col=4)
points(modOrdV, modV, pch=16, col=4)
points(uniqOrd, vel_byOrd_obs, pch=16, col=1)
legend("topleft", c("Observed", "Modeled"), col=c(1,4), pch=16, inset=0.05)



################################################################################
# calculate median discharge by order:
################################################################################
# calculate median discharge (m3/s) by order from hSHEDS flowlines:
uniqOrd = unique(hS$ORDER_)
uniqOrd = sort(uniqOrd)

# remove lakes and non CCC regions from lengths:
hS_noLakes = hS[hS$GLWD_ID == 0 & 
                  hS$gridcode == 0, ]


q_byOrd_obs = rep(NA, length(uniqOrd))
for (i in 1:length(uniqOrd)){
  q_byOrd_obs[i] = median(hS_noLakes$DISCHARGE[hS_noLakes$ORDER_==uniqOrd[i]] )
}

barplot(q_byOrd_obs, names.arg=uniqOrd, 
        main="Outside CCC Median River Discharge",
        xlab="Order",
        ylab="Median River Discharge (cms)")

# extrapolate down to 0-order streams using a POWER LAW:
logQ_lm = lm(log(q_byOrd_obs) ~ uniqOrd)

# determine which orders need to be modeled:
modOrds = c(1:c(min(uniqOrd)-1))
modQ = exp(logQ_lm$coefficients[[2]]*modOrds+logQ_lm$coefficients[[1]])
q_byOrd = c(modS, q_byOrd_obs)

# median  according to the statistical stream order model:
modOrds_all = 1:max(uniqOrd)
modQ_byOrd = exp(logQ_lm$coefficients[[2]]*modOrds_all+logQ_lm$coefficients[[1]])

# plot extrapolation:
plot(c(modOrds, uniqOrd), q_byOrd, type='n',
     main="Outside CCC Median River Discharge by Order",
     xlab="Order",
     ylab="Median Discharge (cms)", log="")
xSeq = seq(-1, max(uniqOrd)+1, length.out=100)
ySeq = exp(logQ_lm$coefficients[[2]]*xSeq+logQ_lm$coefficients[[1]])
lines(xSeq, ySeq, col=4)
points(modOrds, modQ, pch=16, col=4)
points(uniqOrd, q_byOrd_obs, pch=16, col=1)
legend("topleft", c("Observed", "Modeled"), col=c(1,4), pch=16, inset=0.05)



################################################################################
# add attributes to congo hydroSHED shapefile:
################################################################################
# # first use median width and RSSA of each segment based on stream order:
# uniqOrd = unique(hS$ORDER_)
# uniqOrd = sort(uniqOrd)
# 
# hS$width_m = rep(NA, nrow(hS))
# for (i in 1:length(uniqOrd)){
#   ordBoo = hS$ORDER_ == uniqOrd[i]
#   hS$width_m[ordBoo] = modWidth_byOrd[i]
# }
# 
# hS$length_km = hS$LENGTH_KM
# hS$RSSA_km2 = hS$length_km * hS$width_m*1e-3 # convert width from m to km2
# 
# # now add back in the observed median widths:
# mInd = match(hS_OrdGT3$ARCID, hS$ARCID)
# hS$width_m[mInd] = hS_OrdGT3$width_med_
# 
# # zero out lakes and reservoirs:
# hS$width_m[hS$GLWD_ID > 0] = 0
# 
# hS$RSSA_km2 = hS$length_km * hS$width_m*1e-3 # convert width from m to km2
# 
# # rename columns:
# hS$gradient = hS$SLOPE
# hS$velocity_mps = hS$VELOCITY
# hS$discharge_cms = hS$DISCHARGE
# 
# # keep only 6 last columns and create output shapefile: 
# hS = hS[, c((ncol(hS)-5):ncol(hS))]
# 
# write.dbf(hS, hS_oPath)


# create an order table:
ord_tab = as.data.frame(cbind(modOrds_all, 
                              round(width_byOrd, 1), 
                              round(len_byOrd), 
                              round(RSSA_byOrd),
                              round(slope_byOrd, 6),
                              round(vel_byOrd, 3),
                              round(q_byOrd, 4)))
names(ord_tab) = c("Order", "Med_Width_m", "Sum_Len_km", "Sum_Area_km2",
                   "Med_Slope", "Med_Vel_mps", "Med_Q_cms")

write.csv(ord_tab, "E:/research/2018_12_20_Borges_etal_Congo/tabs/ord_tab_congo_CCCexclude.csv", row.names=F)



