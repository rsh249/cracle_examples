library(BIEN)
library(vegdistmod)
library(foreach)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(data.table)
ob.nmin = 5;
t.nmin = 10;
nclus = 16;
clim = stack('sampleclim.gri')
alt = stack('altNA.gri')

if(file.exists('plots.tab')){
	         ext = extent(c(-130, -40, 10, 66));
	#	 plots = read.table('plots.tab', nrows= 13348359); 
		plots = fread('plots.tab', fill=TRUE);
	} else {
	us.plots=BIEN_plot_country('United States', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        can.plots=BIEN_plot_country('Canada', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        mex.plots=BIEN_plot_country('Mexico', cultivated=FALSE, all.taxonomy = TRUE, native.status=TRUE, all.metadata=TRUE)
        ##merge plots
        plots = rbind(us.plots, can.plots,mex.plots);
        write.table(plots, file='plots.tab', append=FALSE, quote=TRUE, sep = "\t", row.names=FALSE)
	plots =	fread('plots.tab', fill=TRUE);
	ext = extent(c(-130, -40, 10, 66))
}



alt.plots = extract(alt, plots[,c('longitude', 'latitude')])
alt.plots = cbind(alt.plots, plots$elevation_m)
altdiff = alt.plots[,1] - alt.plots[,2]
plots = plots[-which(altdiff>500),]


head(plots)
nrow(plots)

#}

###Crop to study area (ext)

	plots = plots[which(plots$lat >= ext[3]),]
	plots = plots[which(plots$lat <= ext[4]),]
	plots = plots[which(plots$lon >= ext[1]),]
	plots = plots[which(plots$lon <= ext[2]),]


#function to quality check lat lon values
decimalplaces = function(x) {
	places=NULL;
	if((x %% 1) != 0){
		places = nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]]);
	} else {
		places = 0;
	}
	return(places);
}


#load('p_works.RData');
##Filter plots table
#plots2 = subset(plots, plots$datasource == 'VegBank'); ##
plots2 = plots;
plots2 = plots2[!is.na(plots2$latitude),];
plots2 = plots2[!is.na(plots2$longitude),];
lat.places = unlist(lapply(plots2$latitude, decimalplaces))

lon.places = unlist(lapply(plots2$longitude, decimalplaces))
plots2 = plots2[which(lat.places>2 & lon.places>2),]
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
		#print(nrow(s));

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
	cl = makeCluster(minclus, type = 'SOCK', outfile = '');
	clusterExport(cl, c("plots2", "ob.nmin"));
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
t_names=unique(plots2$scrubbed_species_binomial)


#Get GBIF data for all taxa
if(!file.exists('./data.tab')){
    extall = getextr(as.character(t_names), clim, schema='flat', factor=8,
    maxrec=5000, rm.outlier=FALSE,
    parallel=TRUE, nclus = nclus,
    nmin = ob.nmin)
    write.table(extall, file='data.tab');
} else {
    extall = read.table('data.tab');
}
head(extall)

extall = na.omit(extall);
ecount = plyr::count(extall$tax)
extall2= extall[which(extall$tax %in% ecount[which(ecount$freq>=ob.nmin),1]),]

extall = extall2;




#PDF construction
cat("Begin PDF construction stage\n");
densall <- dens_obj(extall[,-5],
        clim, manip = 'condi', kern ='gaussian',
        n = 4096, clip = 'range',
        bg.n=200, parallel=TRUE, nclus = nclus);


cat("Calculating names list for all successful taxon models\n");
dens.refnames = NULL;
for (z in 1:length(densall)){
        dens.refnames = c(dens.refnames, as.character(densall[[z]]$name));
}

save.image('test.RData');

predcoll = function(x){
	require(vegdistmod);
	ret = list();
	for(i in 1:length(x)){
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
        	site_ob = cbind('1111', as.character(x[[i]]), lat, lon)
        	site_ob = as.data.frame(site_ob)
       		 colnames(site_ob) = colnames(extall[,1:4])
        	denlist = densall[which(dens.refnames %in% sp)]; 
        	and = and_fun(denlist)
        	opt = get_optim(and);
        	opt.top = apply(opt$origk, 2, mean);
        	actual = raster::extract(clim, cbind(as.numeric(as.character(lon)), as.numeric(as.character(lat))), cellnumbers=TRUE)
        	site_ex = cbind(site_ob, actual)
		ret[[i]] = cbind(site_ex, t(opt.top), length(denlist));
	}
	return(ret);

}




require(parallel);
minclus = nclus/4;
cl = makeCluster(minclus, type="FORK", outfile = '');
#clusterExport(cl, c("plots2", "clim", "extall", "dens.refnames", "densall")); ##Must do if using SOCK cluster type (MPI probably too).
p = proc.time();
#splits = clusterSplit(cl, p_works);
coll2 = parLapply(cl, p_works, predcoll);
proc.time() - p;
stopCluster(cl);



collection = rbind(coll2[[1]][[1]]);
for(k in 2:length(coll2)){
	collection= rbind(collection, coll2[[k]][[1]]);
}








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



collection = collection[which(collection$'length(denlist'>= t.nmin),]

#plot diversity distribution among sites
hist(collection$'length(denlist)')

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
               'bio11', 'MAT', 'MaximumT', 'MinimumT')
for(i in 6:(6+(nlayers(clim)-1))){
  var = collection[,i];  
	var.origkde = collection[,i+nlayers(clim)]

  if(colnames(collection)[i] %in% adjustment){
	#some temperature variables are stored and operated on as Celcius*10. We are reversing that here.

    print(colnames(collection)[i])
    var = var/10;
    var.origkde = var.origkde/10
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

plot(collection[,'length(denlist)'], abs(collection$MAT.origkde - collection$MAT)/10, xlab = 'n taxa', ylab = 'MAT Anomaly Rate (C)')

attach(collection)
plot(MAT, MAT.origkde)
#plot(bio5, bio5.origkde)
abline(0,1)
summary
studvars = c('MAT', 'maxtemp', 'mintemp', 'MAP', 'wbalann', 'WINTERLEN')
MSvars = summary[which(summary[,1]%in% studvars),]
MSvars = MSvars[,-6]
for(i in 2:ncol(MSvars)){
  MSvars[,i] = sprintf("%.3f", as.numeric(MSvars[,i]))
}
colnames(MSvars) = c('var', 'mean error', 'mean absolute anomaly', 'median absolute anomaly', 'RMSE', 'r2', 'pearson', 'spearman')
print(MSvars)
write.table(MSvars, file = 'MSvars.tab', sep='\t')

library(plyr)
df = cbind(abs(collection$MAT - collection$MAT.origkde), collection[,'length(denlist)'])
colnames(df) = c('MAT Median Anomaly (C)','r')
df = as.data.frame(df)
dp=ddply(df, .(cut(df$r, 10)), colwise(mean));
dp[,'MAT Median Anomaly (C)'] = dp[,'MAT Median Anomaly (C)']/10
plot(dp, ylim=c(0,3), xlab = 'N Species', ylab = 'MAT Median Anomaly (C)')

raster::plot(clim[[1]], col = viridis::viridis(100))
points(as.numeric(as.character(collection$lon)), as.numeric(as.character(collection$lat)), pch=20, cex=0.3)

save.image('end.RData');


anomaly.melt = melt(collection, 
	measure.vars=c('MAT', 'maxtemp', 'mintemp', 'MAP', 'wbalann', 'WINTERLEN'),
	variable.name = 'ClimVar',
	value.name = 'Absolute_Anomaly')

ggplot(data=anomaly.melt) +
	geom_histogram(aes(x=Absolute_Anomaly, fill = ClimVar), bins =30) +
	facet_wrap(~ClimVar, scales='free') +
	theme(panel.background=element_rect(fill=NA), legend.title=element_blank())+
	scale_fill_manual(values=brewer.pal('Dark2', n=6))





q('no');



