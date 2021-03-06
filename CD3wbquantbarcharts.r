# Create a bar plot of quantified Western blot band intensities for CD3 subunits using ggplot2.
# TCRzeta, CD3 epsilon, and CD3 gamma are plotted separately.

library(ggplot2)

setwd('C:\\Users\\jmb\\Desktop\\LabFiles\\plcgpaper\\WBs')

# Read in data and make a single data frame including quant data for all three proteins
dat = read.table('CD3wbquantdata.txt', sep='\t', header=TRUE, colClasses = c('numeric','character','character','character','character'))

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

### Make separate plots for each protein
###	TCRz
datTCRz = dat2[which(dat2$prot=='TCRz'),]

datTCRz <- datTCRz[with(datTCRz,order(-geno)), ] ## Sorting
datTCRz$geno <- ordered(datTCRz$geno, levels=levels(datTCRz$geno)[unclass(datTCRz$geno)])

label.df <- data.frame(Group = c('0m', '3m'),
	Value = c(.05, .05))

ggplot(data = datTCRz, aes(x = time, y = ratio, fill = geno, width = 0.6)) +
    geom_bar(aes(fill = geno),group = 'geno', colour='black', stat = 'identity', position=position_dodge(), width=0.7) +
	scale_fill_manual(values=c('white','lightgrey')) +
	xlab('') +
	ylab(expression(paste('Ratio (T cell receptor', zeta,'/ GAPDH)'))) +
	scale_x_discrete(labels = c('0 minutes', '3 minutes')) +
	scale_y_continuous(expand = c(0,0), limits = c(0,65), breaks = seq(0,65,by=5), labels = seq(0,65,by=5)) +
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
###	CD3E
datCD3e = dat2[which(dat2$prot=='CD3E'),]

datCD3e <- datCD3e[with(datCD3e,order(-geno)), ] ## Sorting
datCD3e$geno <- ordered(datCD3e$geno, levels=levels(datCD3e$geno)[unclass(datCD3e$geno)])

label.df <- data.frame(Group = c('0m', '3m'),
	Value = c(.05, .05))

ggplot(data = datCD3e, aes(x = time, y = ratio, fill = geno, width = 0.6)) +
    geom_bar(aes(fill = geno),group = 'geno', colour='black', stat = 'identity', position=position_dodge(), width=0.7) +
	scale_fill_manual(values=c('white','lightgrey')) +
	xlab('') +
	ylab(expression(paste('Ratio (CD3', epsilon,'/ GAPDH)'))) +
	scale_x_discrete(labels = c('0 minutes', '3 minutes')) +
	scale_y_continuous(expand = c(0,0), limits = c(0,65), breaks = seq(0,65,by=5), labels = seq(0,65,by=5)) +
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
	
###	CD3G
datCD3g = dat2[which(dat2$prot=='CD3G'),]

datCD3g <- datCD3g[with(datCD3g,order(-geno)), ] ## Sorting
datCD3g$geno <- ordered(datCD3g$geno, levels=levels(datCD3g$geno)[unclass(datCD3g$geno)])

label.df <- data.frame(Group = c('0m', '3m'),
	Value = c(.05, .05))

ggplot(data = datCD3g, aes(x = time, y = ratio, fill = geno, width = 0.6)) +
    geom_bar(aes(fill = geno),group = 'geno', colour='black', stat = 'identity', position=position_dodge(), width=0.7) +
	scale_fill_manual(values=c('white','lightgrey')) +
	xlab('') +
	ylab(expression(paste('Ratio (CD3', gamma,'/ GAPDH)'))) +
	scale_x_discrete(labels = c('0 minutes', '3 minutes')) +
	scale_y_continuous(expand = c(0,0), limits = c(0,65), breaks = seq(0,65,by=5), labels = seq(0,65,by=5)) +
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
