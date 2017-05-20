# Script to create bar chart of Itk Western blot data

library(ggplot2)

setwd('C:\\Users\\jmb\\Desktop\\LabFiles\\plcgpaper\\WBs')

# read in data and create a data frame. Quantitative data is assumed to be normalized to GAPDH.

dat = read.table('ITKwbquantdata.txt', sep='\t', header=TRUE, colClasses = c('numeric','character','character','character','character'))

std <- function(x) sd(x)/sqrt(length(x))

dat$label = paste(dat$protein, dat$time, dat$geno)

dat2 = as.data.frame(
	do.call(rbind, 
		lapply(unique(dat$label), function(x){
			ratio = mean(dat$ratio[which(dat$label==x)])
			prot = dat$protein[which(dat$label==x)][1]
			time = dat$time[which(dat$label==x)][1]
			geno = dat$geno[which(dat$label==x)][1]
			se = std(dat$ratio[which(dat$label==x)])
			sd = sd(dat$ratio[which(dat$label==x)])
			return(data.frame(ratio,se,sd,prot,time,geno))
		})
	)
)

datITK = dat2[which(dat2$prot=='ITK'),]

datITK <- datITK[with(datITK,order(-geno)), ] ## Sorting
datITK$geno <- ordered(datITK$geno, levels=levels(datITK$geno)[unclass(datITK$geno)])

label.df <- data.frame(Group = c('0m', '3m'),
	Value = c(.05, .05))

ggplot(data = datITK, aes(x = time, y = ratio, fill = geno, width = 0.6)) +
    geom_bar(aes(fill = geno),group = 'geno', colour='black', stat = 'identity', position=position_dodge(), width=0.7) +
	scale_fill_manual(values=c('white','lightgrey')) +
	xlab('') +
	ylab(paste('Ratio (Itk / GAPDH)')) +
	scale_x_discrete(labels = c('0 minutes', '3 minutes')) +
	scale_y_continuous(expand = c(0,0), limits = c(0,.125), breaks = seq(0,.125,by=.025), labels = seq(0,.125,by=.025)) +
	geom_errorbar(aes(ymax = ratio + sd, ymin = ratio - sd), position = position_dodge(width=0.5), width = 0.2) +
#	geom_text(data = label.df, label = '*') +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		#axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
		axis.line.x = element_line(color='black', size=0.6),
		axis.line.y = element_line(color='black', size=0.6),
		axis.text.x = element_text(colour='black', size=19),
		axis.text.y = element_text(colour='black', size=19),
		axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=26),
		axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=26),
		axis.ticks = element_blank(),
		legend.position = c(.85,.94),
		legend.text = element_text(size = 19),
		legend.title=element_blank()
	)
