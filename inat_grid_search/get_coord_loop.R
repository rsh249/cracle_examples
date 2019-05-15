library(raster)
library(cRacle)
library(ggplot2)
library(rinat)
library(rgbif)
nclus=24
ob.nmin=5

#Geographic Bounds lat, lon)
ext = extent(c(-125,-60, 15, 60))
n = 125
#n = 20
#ext = extent(c(-80, -65, 40,50))
ext.larger = extent(ext[1] - 10, ext[2] + 10, ext[3] - 10, ext[4] + 10)

############EDITS REQUIRED
#if you need to get the worldclim data run lines 17:20 FIRST
#getwc2 = 'wget http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_2.5m_bio.zip'
#system(getwc2)
#unz = 'unzip wc2.0_2.5m_bio.zip'
#system(unz)
#####################################

##EDIT PATH TO WORLDCLIM DATA BELOW:$$$$$$$$$$$$$$$$
clim = stack(list.files('/usr/share/data/wc2.0/bio2.5', pattern="*.tif", full.names = T))
#clim = stack(list.files('/usr/share/data/wc2.0/bio0.5', pattern="*.tif", full.names = T))
######################
clim = crop(clim[[1]], ext.larger)
lat = ext[4]
lon = ext[1]

inc.x = (ext[2] - ext[1]) / n
inc.y = (ext[4] - ext[3]) / n
s.lat = vector()
s.lon = vector()
#Generate list of coords
while (lat >= ext[3]) {
  x = c(lat, lon)
  s.lon = append(s.lon, lon)
  s.lat = append(s.lat, lat)
  if (lon == ext[1]) {
    (lon = lon + inc.x)
    (lat = lat - inc.y)
  } else if (lon >= ext[2])  {
    lon = ext[1]
  } else if (lon > ext[1]) {
    lon = lon + inc.x
  }
  print(x)
  
}

points = cbind(s.lon, s.lat)
#plot(clim[[1]])
#points(points, cex=0.2)

#source iqnat function
source('query_function.R')
#loop for points
# library(parallel)
# cl = makeCluster(nclus, type = 'SOCK')
# splits = list()
# for(z in 1:nrow(points)){splits[[i]] = points[i,2:1]}
# par_inat = parSapply(cl, splits, query_inat)
# stopCluster(cl)

# sp_list = list()
# for (i in 1:nrow(points)) {
#   sp_list[[i]] = query_inat(points[i,2:1])
# }

library(parallel)
cl = makeCluster(nclus, type = 'FORK'); #Has to be fork clusters 
p = proc.time();
splits = clusterSplit(cl, 1:nrow(points));
df.splits = list()
for(nn in 1:length(splits)){
  df.splits[[nn]] = points[splits[[nn]], 2:1]
}
sp_list = parSapply(cl, df.splits, query_inat);
proc.time() - p;
stopCluster(cl);



sp_l = sp_list[[1]]
for(s in 2:length(sp_list)){
  sp_l = c(sp_l, sp_list[[s]])
}
sp_list2 = sp_list
sp_list=sp_l
#generate species master list from inat function
sp_list = as.matrix(sp_list)
sp_list = as.data.frame(sp_list)
master = unlist(sp_list$V1)
#remove author names from master list
master = sub('^(\\w+\\s+\\w+).*', '\\1', master)
master = unique(master)

save.image('1.RData')
# Raw and extract after      !!!
if(file.exists('distdata.csv')){extr = read.csv('distdata.csv')} else{
extr = getextr(
  master,
  clim,
  maxrec = 2000,
  schema = 'raw', 
  repo='gbif',
  rm.outlier = FALSE,
  nmin = ob.nmin,
  parallel = T, nclus = nclus
)
write.csv(extr, file='distdata.csv')
}

extr = na.omit(extr);
ecount = plyr::count(extr$tax)
extall2= extr[which(extr$tax %in% ecount[which(ecount$freq>=ob.nmin),1]),]

extr = extall2;





#for each site w/sp data: compute cracle output
densall = dens_obj(
  extr,
  clim,
  manip = 'condi',
  bw = 'nrd0',
  kern = 'optcosine',
  clip = '99conf',
  n = 2048,
  parallel = TRUE,
  nclus = nclus,
  bg.n = 1000
)

cat("Calculating names list for all successful taxon models\n");
dens.refnames = NULL;
for (z in 1:length(densall)){
  dens.refnames = c(dens.refnames, as.character(densall[[z]]$name));
}
points=as.data.frame(points)
site_ex = extract(clim, cbind(points$s.lon, points$s.lat))
ex.points = cbind(points, site_ex)
ex.points$cracle.raw = NA


for (i in 1:nrow(sp_list)) {
#for (i in 1:10) {
   
    sp = sub('^(\\w+\\s+\\w+).*', '\\1', sp_list$V1[[i]])
    sp = unique(sp)
    # denlist is the set of model objects for species in 'sp'
    denlist = densall[which(dens.refnames %in% sp)]; 
    if(length(denlist)<=4){next}
    and = and_fun(denlist)
    opt = get_optim(and);
    opt.top = apply(opt$origk, 2, mean);
    opt.topm = apply(opt$means, 2, mean);
    
    ex.points[i,'cracle.raw'] = opt.top
    ex.points[i,'cracle.means.raw'] = opt.topm
    
    ##MCR/CA code
    sp.names=vector()
    for(qq in 1:length(denlist)){
      sp.names[qq] = as.character(denlist[[qq]]['name']$name)
    }
    subsub = extr[extr$tax %in% sp.names,] #compare sp to sp.names
    mcr.un = MCR(subsub, method='unweight')
    ex.points[i,paste('mcr.un', 'raw', sep='.')] = apply(mcr.un[[1]],2,mean)
    mcr.w = MCR(subsub, method='weight')
    ex.points[i,paste('mcr.w','raw', sep='.')] = apply(mcr.w[[1]],2,mean)
}
#ex.points = na.omit(ex.points)
#saved hits as hits.safe
hits = ex.points[!is.na(ex.points$cracle.raw),]
plot(clim[[1]])
#points(points[,2:1])
points(points, cex = 0.1, pch = 20, col = alpha("black", .2)) #check out cex
points(hits[,1:2], col = 'red', pch =20, cex=0.2)

# hits$highlight = abs(hits$fac4 - hits$site_ex)/10
# 
# hits$highlight[hits$highlight < 1.5] <- NA
# variation = hits
# variation = na.omit(variation)
# points(variation[,1:2], col = 'black', pch =20)



ggplot(data = ex.points, aes(x=abs((ex.points$cracle.raw - ex.points$site_ex)))) +
  geom_histogram() +
  labs(x="Variation from Worldclim, Degrees C", title="Summary Statistics")
ggsave("schema = raw", device = "pdf")

## Loop to test factor of 2,4,6,8,10
sp_sub = sp_list[!is.na(ex.points$cracle.raw),'V1']
t_sub = unlist(sp_list)
extr_sub = extr[extr$tax %in% unique(t_sub), ]
for(fac in c(2,4,6,8,10)){
  subextr = extraction(extr_sub[,1:4], clim=clim, schema = 'flat', factor=fac, rm.outlier=F, nmin=ob.nmin)
  subextr = na.omit(subextr);
  ecount = plyr::count(subextr$tax)
  subextall2= subextr[which(subextr$tax %in% ecount[which(ecount$freq>=ob.nmin),1]),]
  
  subextr = subextall2;
  
  densall = dens_obj(
    subextr,
    clim,
    manip = 'condi',
    bw = 'nrd0',
    kern = 'optcosine',
    clip = '99conf',
    n = 2048,
    parallel = TRUE,
    nclus = nclus,
    bg.n = 1000
  )
  
  cat("Calculating names list for all successful taxon models\n");
  dens.refnames = NULL;
  for (z in 1:length(densall)){
    dens.refnames = c(dens.refnames, as.character(densall[[z]]$name));
  }
  
  
  hits[,paste('cracle', fac, sep='.')] = NA
  
  for (i in 1:length(sp_sub)) {
    print(i)
    #for (i in 1:10) {
    
    sp = sub('^(\\w+\\s+\\w+).*', '\\1', sp_sub[[i]])
    sp = unique(sp)
    sp# denlist is the set of model objects for species in 'sp'
    denlist = densall[dens.refnames %in% sp]; 
    if(length(denlist)<=4){next}
    and = and_fun(denlist)
    opt = get_optim(and);
    opt.kde = apply(opt$origk, 2, mean);
    hits[i,paste('cracle', fac, sep='.')] = opt.kde
    opt.means = opt$means;
    hits[i,paste('cracle.means', fac, sep='.')] = opt.means
    
    
    
    
    
    
    ##MCR/CA code
    sp.names=vector()
    for(qq in 1:length(denlist)){
      sp.names[qq] = as.character(denlist[[qq]]['name']$name)
      }
    subsub = subextr[subextr$tax %in% sp.names,] #compare sp to sp.names
    mcr.un = MCR(subsub, method='unweight')
    hits[i,paste('mcr.un', fac, sep='.')] = apply(mcr.un[[1]],2,mean)
    mcr.w = MCR(subsub, method='weight')
    hits[i,paste('mcr.w', fac, sep='.')] = apply(mcr.w[[1]],2,mean)
  }
  


    
}

ex.points=hits

ex.points$kde.raw = NA
ex.points$kde2 = NA
ex.points$kde4 = NA
ex.points$kde6 = NA
ex.points$kde8 = NA
ex.points$kde10 = NA
ex.points$meansraw = NA
ex.points$means2 = NA
ex.points$means4 = NA
ex.points$means6 = NA
ex.points$means8 = NA
ex.points$means10 = NA
ex.points$wraw = NA
ex.points$w2 = NA
ex.points$w4 = NA
ex.points$w6 = NA
ex.points$w8 = NA
ex.points$w10 = NA
ex.points$uwraw = NA
ex.points$uw2 = NA
ex.points$uw4 = NA
ex.points$uw6 = NA
ex.points$uw8 = NA
ex.points$uw10 = NA

ex.points$kde.raw =  ex.points$cracle.raw - ex.points$site_ex
ex.points$kde2 = ex.points$cracle.2 - ex.points$site_ex
ex.points$kde4 = ex.points$cracle.4 - ex.points$site_ex
ex.points$kde6 = ex.points$cracle.6 - ex.points$site_ex
ex.points$kde8 = ex.points$cracle.8 - ex.points$site_ex 
ex.points$kde10 = ex.points$cracle.10 - ex.points$site_ex
ex.points$meansraw = ex.points$cracle.means.raw - ex.points$site_ex
ex.points$means2 = ex.points$cracle.means.2 - ex.points$site_ex
ex.points$means4 = ex.points$cracle.means.4 - ex.points$site_ex
ex.points$means6 = ex.points$cracle.means.6 - ex.points$site_ex
ex.points$means8 = ex.points$cracle.means.8 - ex.points$site_ex
ex.points$means10 = ex.points$cracle.means.10 - ex.points$site_ex
ex.points$wraw= ex.points$mcr.w.raw - ex.points$site_ex
ex.points$w2 = ex.points$mcr.w.2 - ex.points$site_ex
ex.points$w4 = ex.points$mcr.w.4 - ex.points$site_ex
ex.points$w6 = ex.points$mcr.w.6 - ex.points$site_ex
ex.points$w8 = ex.points$mcr.w.8 - ex.points$site_ex
ex.points$w10 = ex.points$mcr.w.10 - ex.points$site_ex
ex.points$uwraw= ex.points$mcr.un.raw - ex.points$site_ex
ex.points$uw2 = ex.points$mcr.un.2 - ex.points$site_ex
ex.points$uw4 = ex.points$mcr.un.4 - ex.points$site_ex
ex.points$uw6 = ex.points$mcr.un.6 - ex.points$site_ex
ex.points$uw8 = ex.points$mcr.un.8 - ex.points$site_ex
ex.points$uw10 = ex.points$mcr.un.10 - ex.points$site_ex


box_plot <- data.frame(ex.points[1:3], stack(ex.points[28:ncol(ex.points)]))


ggplot(data = box_plot, aes(x= box_plot$ind, (y= abs((box_plot$values))))) +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust = .65))


mea = aggregate(abs(box_plot$values) ~ box_plot$ind, data=box_plot, mean, na.rm=T)
mead = aggregate(abs(box_plot$values) ~ box_plot$ind, data=box_plot, median, na.rm=T)
sdev = aggregate(abs(box_plot$values) ~ box_plot$ind, data=box_plot, sd, na.rm=T)
iqr2 = aggregate(abs(box_plot$values) ~ box_plot$ind, data=box_plot, FUN=quantile, probs=0.25, na.rm=T)
iqr3 = aggregate(abs(box_plot$values) ~ box_plot$ind, data=box_plot, FUN=quantile, probs=0.75, na.rm=T)
summ = cbind(mea, sdev, mead, iqr2, iqr3)
summ = summ[,c(1,2,4,6,8,10)]
colnames(summ) = c('ind', 'mean', 'sdev', 'median', 'iqr2', 'iqr3')

cou = box_plot %>% group_by(ind) %>% count(is.na(values))
summ = cbind(summ, cou[which(cou$`is.na(values)`==FALSE),3])
colnames(summ) = c('ind', 'mean', 'sdev', 'median', 'iqr2', 'iqr3', 'n')

library(ggsci)
test_spdf <- as(clim, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("Degrees", "x", "y")
mapp = ggplot(data=test_df) +  
  geom_tile(aes(x=x, y=y, fill=Degrees)) + 
  coord_equal() + theme_bw() + scale_fill_gsea() +
  labs(x="Longitude", y="Latitude", title='iNaturalist Grid Search') +
  geom_point(data=points, aes(x=s.lon, y=s.lat), alpha=0.01, pch = 15, size=2) +
    geom_point(data=ex.points, aes(x=s.lon, y=s.lat), alpha=1, pch = 15)

mapp 



p1 = ggplot(data = summ, aes(x= ind, y= mean, ymin =  mean-1.96*(sdev/sqrt(n)), ymax =  mean+1.96*(sdev/sqrt(n)))) +
  geom_pointrange() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust = .65)) +
  labs(x="Estimation Type", y="95% Confidence Interval", title='Mean Anomaly Rates')
p2 = ggplot(data = summ, aes(x= ind, y= median, ymin =  iqr2, ymax =  iqr3)) +
  geom_pointrange() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, vjust = .65)) +
  labs(x="Estimation Type", y="Median +/- 1 Quartile", title='Median Anomaly Rates')
p2

p3 = ggplot(data = box_plot, aes(x=abs((box_plot$values)))) +
  geom_histogram(bins=50) + theme_linedraw() +
  labs(x="Variation from Worldclim, Degrees C", y= "Site Count", title="Summary Statistics") +
  facet_wrap(box_plot$ind, nrow = 4) + theme(axis.text.x = element_text()) +
  theme(axis.text.x  = element_text(angle=45, vjust = .65))
p3


gsa = ggarrange(mapp, p1, p2, ncol=1, nrow=3, labels="AUTO")
ggsave(filename='Figure1.png', plot=gsa, device=NULL, width = 4.25, height=9.5, dpi=500)
