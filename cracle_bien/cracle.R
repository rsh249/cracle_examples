
library(raster)
library(parallel)
library(gbm)
library(cRacle)
library(ggplot2)
library(vroom)
library(BIEN)


ob.nmin = 5;
t.nmin = 5;
nclus = 32; #CHECK YOUR HARDWARE>>>ADJUST nclus to be < the total number of CPU cores.

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

if(file.exists('plots.tab')){
	#	plots = data.table::fread('plots.tab', fill=TRUE); 
		plots=vroom::vroom('plots.tab', delim='\t')
	} else {
        plots.TEAM = BIEN_plot_datasource('TEAM', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        plots.VegBank = BIEN_plot_datasource('VegBank', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        plots.CTFS = BIEN_plot_datasource('CTFS', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        plots.gill = BIEN_plot_datasource('gillespie', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        plots.CVS = BIEN_plot_datasource('CVS', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        plots.SALVIAS = BIEN_plot_datasource('SALVIAS', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        plots.FIA = BIEN_plot_datasource('FIA', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        plots = rbind(plots.TEAM, plots.VegBank, plots.CTFS, plots.gill, plots.CVS, plots.SALVIAS, plots.FIA)
        write.table(plots, file='plots.tab', append=FALSE, quote=TRUE, sep = "\t", row.names=FALSE)
        
  	#  plots =	fread('plots.tab', fill=TRUE);
}
ext = extent(c(-160, -10, -57, 80)) #New World


head(plots)
nrow(plots)


decimalplaces = function(x) {
  #function to quality check lat lon values for precision
	places=NULL;
	if((x %% 1) != 0){
		places = nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]]);
	} else {
		places = 0;
	}
	return(places);
}

plots2 = as.data.frame(plots);
plots2 = plots2[!is.na(plots2$latitude),];
plots2 = plots2[!is.na(plots2$longitude),];
#lat.places = unlist(lapply(plots2$latitude, decimalplaces))

#lon.places = unlist(lapply(plots2$longitude, decimalplaces))
#plots2 = plots2[which(lat.places>2 & lon.places>2),]
p_names = unique(plots2$plot_name)
nplots=length(p_names);
optima = NULL;
sitevals=NULL;
lastname='';
#parse plots for suitable candidates:
p_works = list();

filter_spnum2 = function(name){
	names = vector();

	j = 1;
	for(i in 1:length(name)){
		s = plots2[plots2$plot_name==name[[i]],];
		if(length(!is.na(unique(s$scrubbed_species_binomial))) >= ob.nmin 
				& length(unique(s$latitude)) == 1 
					& length(unique(s$longitude)) == 1) { 
			names[[j]] = as.character(name[[i]]);
			j = j+1;	
		}
	}
	return(names)
}


filter_spnum2 = compiler::cmpfun(filter_spnum2);
p_names=unique(plots2$plot_name)

if(file.exists('p_works.txt')){ 
		cat("File p_works.txt exists. OPENING EXISTING VERSION!\n");
		p_works = read.table('p_works.txt')[[1]];
	} else {

	require(parallel);
	minclus = nclus;
	cl = makeCluster(minclus, type = 'FORK', outfile = '');
	#clusterExport(cl, c("plots2", "ob.nmin"));
	#cl = makeCluster(minclus, type ='FORK');
	p = proc.time();
	splits = clusterSplit(cl, p_names);
	p_works2 = parLapply(cl, splits, filter_spnum2);
	proc.time() - p;
	stopCluster(cl);
	p_works = unlist(p_works2);
	write.table(p_works, file='p_works.txt')
}


#get final taxon list:
plots3 = plots2[which(plots2$plot_name %in% p_works),]

plots2=plots3
rm(plots3)

t_names=unique(plots2$scrubbed_species_binomial)

clim = crop(clim, ext)
#Get GBIF data for all taxa
if(!file.exists('./data.tab')){
  extall = cRacle::getextr(
    as.character(t_names),
    clim,
    schema = 'raw',
    maxrec = 5000,
    rm.outlier = FALSE,
    parallel = TRUE,
    nclus = nclus,
    nmin = ob.nmin
  )
    write.table(extall, file='data.tab');
} else {
    extall = read.table('data.tab');
}
head(extall)
save.image('savepoint1.RData')

extall = na.omit(extall);
ecount = plyr::count(extall$tax)
extall2= extall[which(extall$tax %in% ecount[which(ecount$freq>=ob.nmin),1]),]

extall = extall2;

extr = extraction(extall[,1:4], clim, schema='flat', factor = 2, rm.outlier = F, nmin = ob.nmin)
extall = extr; ##

nrow(ecount)
nrow(extall)

#PDF construction
cat("Begin PDF construction stage\n");
densall <- dens_obj(extall[,-5],
        clim, manip = 'condi', kern ='optcosine',
        n = 2048, clip = 'range',
        bg.n=1000, parallel=TRUE, nclus = nclus);


cat("Calculating names list for all successful taxon models\n");
dens.refnames = NULL;
for (z in 1:length(densall)){
        dens.refnames = c(dens.refnames, as.character(densall[[z]]$name));
}

#save.image('test.RData');

predcoll = function(x){
	require(cRacle);
	ret = list();
	for(i in 1:length(x)){
	  out = tryCatch({
		lat = NULL
		lon = NULL
		s = NULL
		sp= NULL
	#	s = subset(plots2, plots2$plot_name==x);
	        s = plots2[plots2$plot_name==as.character(x[[i]]),];

#    	sp = as.character(unique(s$scrubbed_species_binomial))
		sp = unique(s$scrubbed_species_binomial)
	    	lat = unique(s$latitude[!is.na(s$latitude)]);
    		lon = unique(s$longitude[!is.na(s$longitude)]);
    		opt.top = NULL;
    		actual = NULL;
    		site_ex = NULL;
	
    		nsp = sum(!is.na(sp))

		print(x[[i]]);
		print(sp);
		if(length(sp)==0){return(NULL)};
        	site_ob = cbind('1111', as.character(x[[i]]), lat, lon)
        	site_ob = as.data.frame(site_ob)
       		 colnames(site_ob) = colnames(extall[,1:4])
        	denlist = densall[which(dens.refnames %in% sp)]; 
        	and = and_fun(denlist);
        	opt = get_optim(and);
        	opt.top = apply(opt$origk, 2, mean);
        	actual = raster::extract(clim, cbind(as.numeric(as.character(lon)), as.numeric(as.character(lat))), cellnumbers=TRUE)
        	site_ex = cbind(site_ob, actual)
		ret[[i]] = cbind(site_ex, t(opt.top), length(denlist));
	},
		error = function(cond) {
		  # Choose a return value in case of error
		  ret[[i]] = NA
		}, warning = function(cond) {
		  ret[[i]] = NA
		})

	}
	return(ret);

}




require(parallel);
minclus = 12
#cl = makeCluster(minclus, type="SOCK", outfile = '');
#clusterExport(cl, c("plots2", "clim", "extall", "dens.refnames", "densall")); ##Must do if using SOCK cluster type (MPI probably too).
cl = makeCluster(minclus, type = "FORK")
p = proc.time();
splits = clusterSplit(cl, p_works);
coll2 = parLapply(cl, splits, predcoll);
proc.time() - p;
stopCluster(cl);



# collection = rbind(coll2[[1]][[1]]);
# for(k in 2:length(coll2)){
# 	collection= rbind(collection, coll2[[k]][[1]]);
# }

collection = data.table::rbindlist(unlist(coll2, recursive=F))







##FIGURES and SUMMARY STATS::
##Functions

nrmse <- function(v1, v2)
{
  error = abs(v1-v2);
  sqrt(mean(error^2))/max(error)
}
rmse <- function(v1, v2)
{
  error = abs(v1-v2);
  sqrt(mean(error^2))
}
nrmse2 <- function(v1, v2)
{
  error = abs(v1-v2);
  sqrt(mean(error^2))/(max(v1) - min(v1))
}

#Use nmrse2 for commoon NMRSE definition % error relative to empirical range


##backup collection object before conversion and filtering:
collection2 = collection
collection = collection[which(collection[,'length(denlist)'] >= t.nmin),]

collection$lat = as.numeric(as.character(collection$lat));
collection$lon = as.numeric(as.character(collection$lon));

#collection = collection[collection$lon <= -100,] ##Limit to the West ONLY
collection$`length(denlist)` = as.numeric(collection$`length(denlist)`)
#plot diversity distribution among sites
ggplot(collection) + geom_histogram(aes(collection$'length(denlist)'), binwidth = 1)


##Collect lowres data for validation step 2 if desired. NOTE: The aggregation of climate data usually makes CRACLE projections farther from the WorldClim value at 10arcmin for example. (See Harbert and Nixon, 2015 for reverse example -- HiRes climate matches CRACLE better).
#clim10 = aggregate(clim, fact=4, fun = median)
#extr10 = extract(clim10,collection[,4:3])


mea = vector()
mea.ab = vector()
med.ab = vector()
normerr = vector()
normerr2 = vector()
r2 = vector()
pearson = vector()
spearman = vector()
name = vector()

adjustment = c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10',
            'bio11', 'MAT', 'MaximumT', 'MinimumT') #NO LONGER NEEDED WITH WORLDCLIM2.0

collection=as.data.frame(collection)
for(i in 6:(6+(nlayers(clim)-1))){
  var = collection[,i];  
	var.origkde = collection[,i+nlayers(clim)]

  if(colnames(collection)[i] %in% adjustment){
	#some temperature variables are stored and operated on as Celcius*10. We are reversing that here.
	collection[,i]=collection[,i]/10;
	collection[,i+nlayers(clim)] = collection[,i+nlayers(clim)]/10;
    print('adjust')
	  var = collection[,i];
        var.origkde = collection[,i+nlayers(clim)];
  }

  mea[i-5] = mean(var.origkde - var, na.rm=TRUE)
  mea.ab[i-5] = mean(abs(var.origkde - var), na.rm=TRUE)
  med.ab[i-5] = median(abs(var.origkde - var), na.rm=TRUE)
  normerr[i-5] =rmse(var, var.origkde)
  normerr2[i-5] = nrmse2(var, var.origkde)
  r2[i-5] = cor(var, var.origkde) ^ 2
  pearson[i-5] = cor(var, var.origkde)
  spearman[i-5] = cor(var, var.origkde, method='spearman')
  name[i-5] = colnames(collection)[i]

}


summary = cbind(name, mea, mea.ab, med.ab, normerr, normerr2, r2, pearson, spearman)

#replace with ggplot density: plot(collection[,'length(denlist)'], abs(collection$wc2.0_bio_2.5m_01.origkde - collection$wc2.0_bio_2.5m_01), xlab = 'n taxa', ylab = 'MAT Anomaly Rate (C)')
ggplot(collection) + geom_hex(aes(x=collection$`length(denlist)`, y = abs(collection$wc2.0_bio_2.5m_01.origkde - collection$wc2.0_bio_2.5m_01)))
attach(collection)
plot(wc2.0_bio_2.5m_01, wc2.0_bio_2.5m_01.origkde, cex=0.1)
#plot(bio5, bio5.origkde)
abline(0,1)
summary

plot(clim[[1]], col=viridis::magma(999))
points(collection$lon, collection$lat, pch = 20, cex =0.2)

#collsafe = collection
anomaly = cbind(collection[,1:5], 
                abs(collection[,6:(nlayers(clim)+5)] - collection[,(nlayers(clim)+6):(2*nlayers(clim)+5)]), 
                collection[,ncol(collection)])
library(reshape2)
anomaly.melt = melt(anomaly, 
	#measure.vars=c('MAT', 'maxtemp', 'mintemp', 'MAP', 'wbalann', 'WINTERLEN'),
	measure.vars=names(clim)[c(1,5,6,12,15,19)],
	variable.name = 'ClimVar',
	value.name = 'Absolute_Anomaly')

save.image('end.RData');

# Correction models
y = collection; #data
sam = sample(nrow(y)/2, replace = F)
train.col = y[sam,]
train.col = dplyr::sample_frac(y, 0.5, replace = FALSE)
test.col = y[-which(y$tax %in% train.col$tax),]



# ##SVR Model testing and prediction
# best.mod = list()
# best.err = list()
# 
# for(var in 1:19) {
#   library(foreach)
#   library(doParallel)
#   library(e1071)
#   
#   x.r = names(clim[[var]])
#   x.p = paste(x.r, 'origkde', sep='.')
#   x = as.formula(paste(x.r, " ~ ", x.p)); #model
#   print(x)
#   y = train.col; #data
#   t = test.col
# 
#   gamma = 2 ^ (-1:4)
#   cost = 2 ^ (0:6)
#   epsilon = seq(0.1, 1, 0.3)
#   gam = list()
#   cos = list()
#   eps = list()
#   t.n = 1
#   for (n in 1:length(gamma)) {
#     for (o in 1:length(cost)) {
#       for (q in 1:length(epsilon)) {
#         #tuning[[t.n]] = paste('gamma=', gamma[n], ", cost=", cost[o])
#         gam[[t.n]] = gamma[n]
#         cos[[t.n]] = cost[o]
#         eps[[t.n]] = epsilon[q]
#         t.n = t.n + 1
#       }
#     }
#   }
#   #print(gam)
#   sampleCount = length(gam)
# 
#   registerDoParallel(cores = 32)
#   
#   modelDataSvm <- foreach(i = 1:sampleCount) %dopar% {
#     library(e1071)
#     svm(
#       x,
#       y,
#       probability = TRUE,
#       cross = 0,
#       gamma = gam[[i]],
#       cost = cos[[i]]
#     )
#   }
#   
#   best.pred = lapply(modelDataSvm, predict, t)
#   
#   med.err.fit = vector()
#   for (i in 1:sampleCount) {
#     med.err.fit[i] = median(abs(best.pred[[i]] - t[,x.r]))
#     
#   }
#   
#   plot(med.err.fit, cos)
#   best = which(med.err.fit == min(med.err.fit))
#    if(length(best) > 1){best = best[[1]]}
#   gam[best]
#   cos[best]
#   
#   #Evaluate each tuning for lowest testing error
#   
#   best.pred = predict(modelDataSvm[[best]], t)
#   best.mod[[var]] = modelDataSvm[[best]]
#   best.err[[var]] = med.err.fit[best]
#   plot(best.pred,
#        t[,x.r],
#        pch = 20,
#        cex = 0.5)
#   points(
#     t[,x.p],
#     t[,x.r],
#     col = 'red',
#     cex = 0.5
#   )
#   abline(0, 1, lwd = 2, col = 'purple')
#   
# }

### GBMs
gbms = list()
gbm.mea.err = list()
gbm.med.err = list()
gbm.mea.abs.err = list()
gbm.med.abs.err = list()
gbm.normerr = list()
gbm.normerr2 = list()
gbm.r2 = list()
gbm.pearson = list()
gbm.spearman = list()
gbm.coll = list()
gbm.test = list()
zbest.iter = list()

colnames(train.col)[ncol(train.col)] = 'ntax'
colnames(test.col)[ncol(train.col)] = 'ntax'

for(var in 1:19) {
  library(gbm)
  
  x.r = names(clim[[var]])
  x.p = paste(x.r, 'origkde', sep='.')
  x.n = 'ntax'
  x = as.formula(paste(x.r, " ~ ", x.p, ' + ', x.n)); #model
  print(x)
  y = train.col; #data
  t = test.col
  z = y[,c(x.r, x.p, x.n)]
  gbm.test = gbm(
    x,
    data = z,
    n.trees = 5000,
    distribution = "laplace",
    shrinkage = 0.1,
    interaction.depth = 3,
    bag.fraction = 0.5,
    train.fraction = 0.5,
    n.minobsinnode = 10,
    keep.data = TRUE,
    cv.folds=5, n.cores=5,
    verbose = TRUE
  )
  zbest.iter[[var]] = gbm.perf(gbm.test, plot.it=F, method='cv')
  gbm.pred = predict(gbm.test, t, n.trees = zbest.iter[[var]][1], type = "link")  
  plot(t[,x.r], gbm.pred,
       pch = 20, col = alpha('black', 0.2),
       cex = 0.5)
  points(
    t[,x.r],
    t[,x.p],
    col = alpha('red', 0.2),
    cex = 0.5
  )
  abline(0,1)
  print(median(abs(t[,x.r] - t[,x.p])))
  print(median(abs(t[,x.r] - gbm.pred)))
  print(cor.test(t[,x.r], t[,x.p]))
  print(cor.test(t[,x.r], gbm.pred))
  
  gbms[[var]] = gbm.test
  gbm.mea.err[[var]] = mean(t[,x.r] - gbm.pred)
  gbm.med.err[[var]] = median(t[,x.r] - gbm.pred)
  gbm.mea.abs.err[[var]] = mean(abs(t[,x.r] - gbm.pred))
  gbm.med.abs.err[[var]] = median(abs(t[,x.r] - gbm.pred))
  gbm.normerr[[var]] = rmse(t[,x.r], gbm.pred)
  gbm.normerr2[[var]] = nrmse2(t[,x.r], gbm.pred)  
  gbm.r2[[var]] = cor(t[,x.r], gbm.pred)^2
  gbm.pearson[[var]] = cor(t[,x.r], gbm.pred)
  gbm.spearman[[var]] = cor(t[,x.r], gbm.pred, method='spearman')
  gbm.coll[[var]] = gbm.pred
  gbm.test[[var]] = t[,x.p]
}


summary.gbm = cbind(
  summary[, 'name'],
  unlist(gbm.mea.err),
  unlist(gbm.mea.abs.err),
  unlist(gbm.med.abs.err),
  unlist(gbm.normerr),
  unlist(gbm.normerr2),
  unlist(gbm.r2),
  unlist(gbm.pearson),
  unlist(gbm.spearman)
)

colnames(summary.gbm) = colnames(summary)

write.csv(summary.gbm, 'summary.gbm')
write.csv(summary, 'summary')
coll.gbm = gbm.coll[[1]]
for(z in 2:19){
  coll.gbm = cbind(coll.gbm, gbm.coll[[z]])
}

colnames(coll.gbm) = paste(names(clim), ".gbm", sep='')

coll.gbm = cbind(test.col[,1:24], coll.gbm, test.col[,ncol(test.col)])

anoma.gbm = cbind(coll.gbm[,1:5], 
                (coll.gbm[,6:(nlayers(clim)+5)] - coll.gbm[,(nlayers(clim)+6):(2*nlayers(clim)+5)]), 
                coll.gbm[,ncol(coll.gbm)])

# coll.test = gbm.test[[1]]
# for(z in 2:19){
#   coll.test = cbind(coll.test, gbm.test[[z]])
# }
coll.test = test.col
#colnames(coll.test) = paste(names(clim), ".test", sep='')

coll.test = cbind(test.col[,1:24], coll.test, test.col[,ncol(test.col)])

anoma.test = cbind(test.col[,1:5], 
                  (coll.test[,6:(nlayers(clim)+5)] - test.col[,(nlayers(clim)+6):(2*nlayers(clim)+5)]), 
                  test.col[,ncol(test.col)])

library(reshape2)
anomaly.melt.gbm = melt(anoma.gbm, 
                    measure.vars=names(clim)[c(1,5,6,12,17,18)],
                    variable.name = 'ClimVar',
                    value.name = 'Absolute_Anomaly')
anomaly.melt.test = melt(anoma.test, 
                        measure.vars=names(clim)[c(1,5,6,12,17,18)],
                        variable.name = 'ClimVar',
                        value.name = 'Absolute_Anomaly')
anomaly.melt.test = cbind(anomaly.melt.test , rep('CRACLE', nrow(anomaly.melt.test )))
anomaly.melt.gbm = cbind(anomaly.melt.gbm, rep('CRACLE + GBR', nrow(anomaly.melt.gbm)))
colnames(anomaly.melt.test )[ncol(anomaly.melt.test )] = 'model'
colnames(anomaly.melt.gbm)[ncol(anomaly.melt.gbm)] = 'model'
colnames(anomaly.melt.gbm)[19] = 'n.tax'
colnames(anomaly.melt.test )[19] = 'n.tax'
anomaly.merge = rbind(anomaly.melt.test, anomaly.melt.gbm)

library(RColorBrewer)
library(ggsci)
err1 = ggplot(data=anomaly.merge) +
  geom_density(aes(x=Absolute_Anomaly, col = model)) +
  facet_wrap(~ClimVar, scales='free',ncol=1) +  theme_linedraw() +
  #theme(panel.background=element_rect(fill=NA), legend.title=element_blank())+
 # scale_fill_manual(values=brewer.pal('Dark2', n=6)) +
  scale_color_npg() + theme(legend.position="bottom", legend.direction="vertical") +
  ylab('Proportion of CRACLE results') +
  xlab('Model Estimate - WorldClim')
err1
ggsave('cracle_anomalies_bien.png', plot=err1, width=7.25, height=3, units='in', dpi=600)


##scatterplots for top vars showing raw cracle and corrected estimates vs worldclim
head(coll.gbm)
head(collection)

nl = nlayers(clim)

coll.catch = data.frame(stringsAsFactors = FALSE)

gbm.catch = data.frame(stringsAsFactors = FALSE)


last = 1
measure.vars = c(1,5,6,12,17,18)
for(nn in measure.vars){
  nam = names(clim[[nn]])
  midp = cbind(test.col[,1:5], test.col[,nam], test.col[,paste(nam, 'origkde', sep ='.')], rep(nam, nrow(test.col)), rep('CRACLE', nrow(test.col)))
  coll.catch= rbind(coll.catch, midp)
  midp.g = cbind(coll.gbm[,1:5], coll.gbm[,nam], coll.gbm[,paste(nam, 'gbm', sep ='.')], rep(nam, nrow(coll.gbm)), rep('CRACLE+GBM', nrow(coll.gbm)))
  gbm.catch= rbind(gbm.catch, midp.g)


  
}


colnames(coll.catch) = c('ind_id', 'tax', 'lat', 'lon', 'cells', 'WorldClim', 'Estimate', 'Parameter', 'Model')
colnames(gbm.catch) = c('ind_id', 'tax', 'lat', 'lon', 'cells', 'WorldClim', 'Estimate', 'Parameter', 'Model')

p.cracle = ggplot(data=coll.catch) + 
  geom_point(aes(x=WorldClim, y =Estimate, colour=factor(Parameter)), size =0.6, alpha=0.7) + 
  facet_wrap(~Parameter, scales = 'free', ncol=1) +
  scale_color_npg() +
  theme_linedraw() + ylab('CRACLE Estimate') + 
  geom_abline(intercept = 0, slope = 1) +
  theme(legend.position='none') + theme(axis.text.x  = element_text(angle=45, vjust = .65))

gbr.cracle = ggplot(data=gbm.catch) + 
  geom_point(aes(x=WorldClim, y =Estimate, colour=factor(Parameter)), size =0.6, alpha=0.7) + 
  facet_wrap(~Parameter, scales = 'free', ncol=1) +
  scale_color_npg() + ylab('CRACLE+GBR Estimate') + 
  geom_abline(intercept = 0, slope = 1) +
  theme_linedraw() + theme(legend.position='none') + theme(axis.text.x  = element_text(angle=45, vjust = .65))


library(ggsci)
test_spdf <- as(clim[[1]], "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("Degrees", "x", "y")
mapp = ggplot(data=test_df) +  
  geom_tile(aes(x=x, y=y, fill=Degrees)) + 
  coord_equal() + theme_bw() + scale_fill_gsea() +
  labs(x="Longitude", y="Latitude", title='BIEN Vegetation Plots') +
  geom_point(data=collection, aes(x=lon, y=lat), alpha=1, pch = 20, size=0.4) 

mapp 


#library(ggpubr)
#sa = ggarrange(p.cracle, gbr.cracle, err1, ncol=3, nrow=1, labels="AUTO")


library("gridExtra")
library("cowplot")
# Arrange plots using arrangeGrob
# returns a gtable (gt)
gt <- arrangeGrob(mapp, p.cracle, gbr.cracle, err1,
                  ncol = 3, nrow = 2, heights=c(1,3),
                  layout_matrix = rbind(c(1,1,1), c(2,3,4)))

p <- ggpubr::as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0.36, 0.69), y = c(1, 0.75, 0.75, 0.75)) 

ggsave(filename='Figure2.png', plot=p, device=NULL, width = 7.25, heigh=9.5, dpi=500)

