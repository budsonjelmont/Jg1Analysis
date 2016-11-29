library(ggplot2)

setwd('C:/Users/Judson/Documents/plcgpaper/Lck substrates for clustering/uncentered correlation avg')

dat = read.xls('PLCr_rebuild_1FDRpTyr_081216 All Lck targets only PPsite and Scansite extended with all protein sites-V2.xls')

#remove asterisks from the phosphosites
psiteCol = which(colnames(dat)=='phosphosite.annotated')
psites = lapply(dat[,psiteCol],function(r){
			gsub('\\*','',r)
		})
dat[,psiteCol] = unlist(psites)

#log2 transform
start = which(colnames(dat) == 'SILAC.ratio.32.for.user.selected.SILAC.timepoint1')
end = which(colnames(dat) == 'SILAC.ratio.32.for.user.selected.SILAC.timepoint6')
dat[,start:end] = log2(dat[,start:end])

#scale data
dat[,start:end] = scale(dat[,start:end])

#make ids_phosphosite
idCol = which(colnames(dat)=='timecourse.manual.naming..protein.name.manual')
dat$id_psite = apply(dat,1,function(x){
	paste(x[idCol],x[psiteCol],sep='_')
})



NAcount = apply(dat, 1, function(x){
	sum(is.na(x[start:end]))}
 )

dat$NAcount = NAcount

dat2 = dat[-which(dat$NAcount > 1),]

#Sort by name
dat2 = dat2[order(dat2$timecourse.manual.naming..protein.name.manual, decreasing = TRUE),]

colnames(dat2)[start:end] = c('0 min', '1 min', '2 min', '3 min', '5 min', '10 min')

#define clusters
# c1 = c('Epsin 2_S172Y186', 'Epsin 2_S173Y186', 'Protein tyrosine phosphatase, non-receptor type 6_S556S557Y564','Glucocorticoid receptor DNA binding factor 1_Y1105',
  # 'CD5_S447Y453','Protein tyrosine phosphatase, non-receptor type 6_S557Y564', 'Protein tyrosine phosphatase, non-receptor type 6_S556Y564',
  # 'ZAP70_Y319')
  
# c2 = c('PAG_S354Y359', 'PAG_S353Y359', 'Lck_Y505', 'Protein tyrosine phosphatase, non-receptor type 6_Y64',
  # 'Pyruvate dehydrogenase E1 alpha subunit, testis specific form_Y299', 'Cold inducible RNA binding protein_Y164', 'Cytoskeleton associated protein 1_Y107',
  # 'HPK1_Y28', 'Epsin 2_S192S195Y196', 'Epsin 2_S192Y196T200', 'Epsin 2_Y186', 'Epsin 2_Y196','Epsin 2_S195Y196S199', 'PAG_Y227S229',
  # 'Pyruvate dehydrogenase E1 alpha subunit, testis specific form_Y287Y299', 'Tyrosine-protein kinase SgK223_Y411',
  # 'Pyruvate dehydrogenase E1 alpha subunit, testis specific form_S293Y299', 'Pyruvate dehydrogenase E1 alpha subunit, testis specific form_Y287S298',
  # 'Pyruvate dehydrogenase E1 alpha subunit, testis specific form_S291Y299', 'T ciell antigen receptor, zeta_Y123', 'T cell antigen receptor, zeta_Y83',
  # 'CD5_T445Y453', 'T cell antigen receptor, zeta_Y142S146', 'T cell antigen receptor, zeta_Y142T147', 'Protein tyrosine phosphatase, non-receptor type 6_Y564',
  # 'Spectrin, alpha, non-erythrocytic 1 (alpha-fodrin)_Y1579', 'CD5_Y453', 'ZAP70_Y493', 'ZAP70_Y292', 'PAG_Y227', 'T cell antigen receptor, zeta_Y142',
  # 'T cell antigen receptor, zeta_Y111', 'ZAP70_S289Y292', 'ZAP70_T286Y292', 'ZAP70_Y492', 'PAG_Y359', 'ZAP70_Y492Y493', 'ZAP70_S491Y493', 'ZAP70_S491Y492',
  # 'CD5_S439Y453', 'Protein tyrosine phosphatase, non-receptor type 6_T555Y564', 'Potassium channel, voltage gated, shaker related subfamily beta, member 2_T18Y25',
  # 'Potassium channel, voltage gated, shaker related subfamily beta, member 2_S20Y25', 'Dynamin binding protein_Y515', 'DEAD (Asp-Glu-Ala-Asp) box polypeptide 17_Y580',
  # 'ZAP70_Y598', 'Protein tyrosine phosphatase, non-receptor type 6_Y98', 'PAG_Y417', 'ZAP70_T300Y319')
# c3 = c('T cell antigen receptor, zeta_Y153', 'Cold inducible RNA binding protein_Y167', 'ZAP70_Y493T494', 'ZAP70_Y492T494', 'Lipopolysaccharide specific response protein 5_Y100',
  # 'FYB_S558Y559','Epsin 2_S192Y196')
# c4 = c('Chromosome 11 open reading frame 59_Y138', 'Chromosome 11 open reading frame 59_Y140', 'Gamma synergin_Y745', 'HIP55_Y140', 'Junction mediating and regulatory protein_Y133',
  # 'ZAP70_Y597', 'VAV1_Y541', 'Talin_Y71', 'Lck_Y192', 'T cell antigen receptor, zeta_Y72', 'Cyclin dependent kinase 2_Y19', 'UPF3B_Y442', 'CDC2_Y15Y19', 'Cytoskeleton associated protein 1_Y114',
  # 'ZAP70_Y597Y598', 'Serine/threonine protein kinase Nek7_Y201', 'Lck_Y470', 'Cold inducible RNA binding protein_Y142', 'HIP55_Y224', 'VAV1_Y791', 'RNA binding protein 4_Y345')

setwd('C:/Users/Judson/Documents/plcgpaper/Lck substrates for clustering/uncentered correlation avg/k6 clusters')

clust = read.table('cluster 6.txt', sep='\t', header=TRUE)
clust = clust[-1,]
avgProfile = unlist(lapply(clust[,2:ncol(clust)], mean, na.rm=TRUE))
clust$id = as.character(clust$id)
clust$datType = 0
clust$lwid = 0.001
clust = rbind(clust, c('avgProfile', avgProfile, 1,0.75))

###make data frame
nr = nrow(clust)
lastcol = ncol(clust)-2
dat = data.frame(timePoint = as.numeric(as.character(c(rep(0,nr),rep(1,nr),rep(2,nr),rep(3,nr),rep(5,nr),rep(10,nr)))),
  # peptide = c(rep(dat$id
	# lapply(dat$id, rep, nc-1),'avg'),
  peptide = rep(clust$id,6),
  ratio = as.numeric(as.character(unlist(clust[,2:lastcol]))),
  type = rep(unlist(clust$datType),6),
  lwid = rep(unlist(clust$lwid),6)
  )

#factor peptide ids by an integer specifying the order elements should be plotted
#to ensure that avg profile line is placed on top
dat$group = as.factor(apply(format(dat[,c('type', 'peptide')]), 1, paste, collapse=' '))

####

r = range(na.omit(dat$ratio))
r[1] = floor(r)[1] 
r[2] = ceiling(r)[2] 

ggplot(dat, aes(x=timePoint, y=ratio)) +
  xlab('\nTime (minutes)') +
  ylab('\nlog2(Jgamma1/Jgamma.WT)') +
#  geom_line(aes(y=ratio, group=group, color=type), size = lwid) +
  geom_line(aes(y=ratio, group=group, color=type, size = lwid)) +
  scale_color_manual(values=c('lightgray','red')) +
  scale_x_continuous(breaks=c(0,5,10)) +
  scale_y_continuous(breaks = seq(r[1],r[2],by=1), labels = seq(r[1],r[2],by=1)) +
  geom_hline(aes(yintercept=0), linetype='dashed', size=0.75) +
#  scale_size(range=c(0.75,1.5), guide=FALSE) +	
#  ggtitle('Zap70 Y315') +
  theme(legend.title = element_blank(),
    legend.key = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    #axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
    axis.line.x = element_line(color='black', size=0.6),
    axis.line.y = element_line(color='black', size=0.6),
    axis.text.x = element_text(colour='black', size=19),
    axis.text.y = element_text(colour='black', size=19),
    axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=26),  #I don't think this is doing anything
    axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=26),  #ditto
    legend.position = 'none')

