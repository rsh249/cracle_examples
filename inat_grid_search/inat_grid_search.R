library(raster)
library(cRacle)
library(ggplot2)
library(rinat)
library(rgbif)
nclus=30
ob.nmin=5

#Geographic Bounds lat, lon)
ext = extent(c(-100,-60, 15, 60))
ext.larger = extent(ext[1] - 10, ext[2] + 10, ext[3] - 10, ext[4] + 10)
clim = stack(list.files('/usr/share/data/wc2.0/bio2.5', pattern="*.tif", full.names = T))
#clim = stack(list.files('/usr/share/data/wc2.0/bio0.5', pattern="*.tif", full.names = T))

clim = crop(clim[[1]], ext.larger)
lat = ext[4]
lon = ext[1]
n = 100
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
plot(clim[[1]])
points(points)

#source inat function
source('~/cRackle testing/query_function.R')
#loop for points
# library(parallel)
# cl = makeCluster(nclus, type = 'SOCK')
# splits = list()
# for(z in 1:nrow(points)){splits[[i]] = points[i,2:1]}
# par_inat = parSapply(cl, splits, query_inat)
# stopCluster(cl)

sp_list = list()
for (i in 1:nrow(points)) {
  sp_list[[i]] = query_inat(points[i,2:1])
}

#generate species master list from inat function
sp_list = as.matrix(sp_list)
sp_list = as.data.frame(sp_list)
master = unlist(sp_list$V1)
#remove author names from master list
master = sub('^(\\w+\\s+\\w+).*', '\\1', master)
master = unique(master)

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



#collect master list of species
#query gbif for all

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
    ex.points[i,'cracle.raw'] = opt.top
    
}
#ex.points = na.omit(ex.points)

#Statistics



hits = na.omit(ex.points)
plot(clim[[1]])
points(points[,2:1])
points(points) #check out cex
points(hits[,1:2], col = 'red', pch =20)
hits$highlight = abs(hits$fac4 - hits$site_ex)/10

hits$highlight[hits$highlight < 1.5] <- NA
variation = hits
variation = na.omit(variation)
points(variation[,1:2], col = 'black', pch =20)




ggplot(data = ex.points, aes(x=abs((ex.points$cracle.raw - ex.points$site_ex)/10))) +
  geom_histogram() +
  labs(x="Variation from Worldclim, Degrees C", title="Summary Statistics")
ggsave("schema = raw", device = "pdf")

## Loop to test factor of 2,4,6,8,10

for(fac in c(2,4,6,8,10)){
  subextr = extraction(extr[,1:5], clim=clim, schema = 'flat', factor=fac, rm.outlier=F, nmin=ob.nmin)
  subextr = na.omit(subextr);
  ecount = plyr::count(subextr$tax)
  subextall2= extr[which(subextr$tax %in% ecount[which(ecount$freq>=ob.nmin),1]),]
  
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
  
  
  ex.points[,paste('cracle', fac, sep='.')] = NA
  
  for (i in 1:nrow(sp_list)) {
    print(i)
    #for (i in 1:10) {
    
    sp = sub('^(\\w+\\s+\\w+).*', '\\1', sp_list$V1[[i]])
    sp = unique(sp)
    sp# denlist is the set of model objects for species in 'sp'
    denlist = densall[which(dens.refnames %in% sp)]; 
    if(length(denlist)<=4){next}
    and = and_fun(denlist)
    opt = get_optim(and);
    opt.kde = apply(opt$origk, 2, mean);
    ex.points[i,paste('cracle', fac, sep='.')] = opt.kde
    opt.means = opt$means;
    ex.points[i,paste('cracle.means', fac, sep='.')] = opt.means
    
    ##MCR/CA code
    sp.names=vector()
    for(qq in 1:length(denlist)){sp.names[qq] = as.character(denlist[[qq]]['name']$name)}
    subsub = subset(subextr, subextr$tax %in% sp.names) #compare sp to sp.names
    mcr.un = MCR(subsub, method='unweight')
    ex.points[i,paste('mcr.un', fac, sep='.')] = apply(mcr.un[[1]],2,mean)
    mcr.w = MCR(subsub, method='weight')
    ex.points[i,paste('mcr.w', fac, sep='.')] = apply(mcr.w[[1]],2,mean)
  }
  


    
}


#Create boxplot table. Make it a matrix to begin with. We know the number of col (2) and we know the number of rows (coords * # of col)
#Fill the col with the data and the name associated with the data
#
# numfactor=6
# numrows = nrow(ex.points)*numfactor
# err.df = data.frame(rep(NA, numrows), rep(NA, numrows))
# end=0
# start=0
# for(z in 4:(4+numfactor)){
#   n = nrow(ex.points)
#   start=(1*(z-3))+(n*(z-4))
#   #start=(1*(z-3))+end
#   end=start+n
#   err.df[,1] = ex.points[,3] - ex.points[,z]
#   err.df[,2] = rep(colnames(ex.points)[z], n)
# }

ex.points$raw = NA
ex.points$fac2 = NA
ex.points$fac4 = NA
ex.points$fac6 = NA
ex.points$fac8 = NA
ex.points$fac10 = NA

ex.points$fac2 = ex.points$cracle.2 - ex.points$site_ex
ex.points$fac4 = ex.points$cracle.4 - ex.points$site_ex
ex.points$fac6 = ex.points$cracle.6 - ex.points$site_ex
ex.points$fac8 = ex.points$cracle.8 - ex.points$site_ex 
ex.points$fac10 = ex.points$cracle.10 - ex.points$site_ex
ex.points$raw =  ex.points$cracle.raw - ex.points$site_ex

box_plot <- data.frame(ex.points[1:3], stack(ex.points[26:ncol(ex.points)]))

ggplot(data = box_plot, aes(x= box_plot$ind, (y= abs((box_plot$values)/10)))) +
  geom_boxplot()


ggplot(data = box_plot, aes(x=abs((box_plot$values / 10)))) +
  geom_histogram() +
  labs(x="Variation from Worldclim, Degrees C", title="Summary Statistics") +
  facet_grid(box_plot$ind) + theme(axis.text.x = element_text())

ggsave("factors", device = "pdf")






#colnames(site_ob) = colnames(extall[,1:4])

# sp is list of species from one location

# denlist is the set of model objects for species in 'sp'
# denlist = densall[which(dens.refnames %in% sp)]; 
# and = and_fun(denlist)
# opt = get_optim(and);
# opt.top = apply(opt$origk, 2, mean);

# actual is the WorldClim estimated climate from a location
# actual = raster::extract(clim, cbind(as.numeric(as.character(lon)), as.numeric(as.character(lat))), cellnumbers=TRUE)
# site_ex = cbind(site_ob, actual)

#return an object that includes site_ex, opt.top, and # of species
# ret[[i]] = cbind(site_ex, t(opt.top), length(denlist));

#collect performance stats


#add col for world clim and crackle 


