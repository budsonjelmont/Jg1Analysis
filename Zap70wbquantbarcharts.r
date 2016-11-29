library(ggplot2)

setwd('C:\\Users\\jmb\\Desktop\\LabFiles\\plcgpaper\\WBs')

dat = read.table('Zap70wbquantdata.txt', sep='\t', header=TRUE, colClasses = c('numeric','character'))

std <- function(x) sd(x)/sqrt(length(x))

dat2 = data.frame(avgratio = c(
	  mean(dat[1:3,1]), mean(dat[4:6,1]), mean(dat[7:9,1]),
	  mean(dat[10:12,1]), mean(dat[13:15,1]), mean(dat[16:18,1]),
	  mean(dat[19:21,1]), mean(dat[22:24,1])), 
	 se = c(
	  std(dat[1:3,1]), std(dat[4:6,1]), std(dat[7:9,1]),
	  std(dat[10:12,1]), std(dat[13:15,1]), std(dat[16:18,1]),
	  std(dat[19:21,1]), std(dat[22:24,1])),
	 geno = rep(c('Jgamma1.WT','Jgamma1.WT','Jgamma1','Jgamma1'),2),
	 time = rep(c('0m','3m','0m','3m'),2),
	 site = c(rep('pY319',4), rep('pY493',4))
	)

dat2$label<-factor(raw[,2],levels=c("0m","3m"),ordered=FALSE)

#dat2$label = paste(dat2$geno, dat2$site,sep='_')


###separate plots for each residue
###	pY319
dat319 = dat2[1:4,]

dat319 <- dat319[with(dat319,order(-geno)), ] ## Sorting
dat319$geno <- ordered(dat319$geno, levels=levels(dat319$geno)[unclass(dat319$geno)])

label.df <- data.frame(Group = c('0m', '3m'),
	Value = c(.05, .05))

ggplot(data = dat319, aes(x = time, y = avgratio, fill = geno, width = 0.6)) +
    geom_bar(aes(fill = geno),group = 'geno', colour='black', stat = 'identity', position=position_dodge(), width=0.7) +
	scale_fill_manual(values=c('white','lightgrey')) +
	xlab('Zap70 pY319') +
	ylab('Ratio Zap70 pY319/total Zap70') +
	scale_x_discrete(labels = c('0 minutes', '3 minutes')) +
	scale_y_continuous(expand = c(0,0), limits = c(0,0.35), breaks = seq(0,0.35,by=0.05), labels = seq(0,0.35,by=0.05)) +
	geom_errorbar(aes(ymax = avgratio + se, ymin = avgratio - se), position = position_dodge(width=0.5), width = 0.2) +
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
		axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=26),  #I don't think this is doing anything
		axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=26),  #ditto
		axis.ticks = element_blank(),
		legend.position = c(.25,.9),
		legend.text = element_text(size = 19),
		legend.title=element_blank()
	)

###	pY493
dat493 = dat2[5:8,]

dat493 <- dat493[with(dat493,order(-geno)), ] ## Sorting
dat493$geno <- ordered(dat493$geno, levels=levels(dat493$geno)[unclass(dat493$geno)])

label.df <- data.frame(Group = c('0m', '3m'),
	Value = c(.05, .05))

ggplot(data = dat493, aes(x = time, y = avgratio, fill = geno, width = 0.6)) +
    geom_bar(aes(fill = geno),group = 'geno', colour='black', stat = 'identity', position=position_dodge(), width=0.7) +
	scale_fill_manual(values=c('white','lightgrey')) +
	xlab('Zap70 pY493') +
	ylab('Ratio Zap70 pY493/total Zap70') +
	scale_x_discrete(labels = c('0 minutes', '3 minutes')) +
	scale_y_continuous(expand = c(0,0), limits = c(0,0.35), breaks = seq(0,0.35,by=0.05), labels = seq(0,0.35,by=0.05)) +
	geom_errorbar(aes(ymax = avgratio + se, ymin = avgratio - se), position = position_dodge(width=0.5), width = 0.2) +
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
		axis.title.x = element_text(margin=margin(t=0, r=0, b=10.8, l=0), size=26),  #I don't think this is doing anything
		axis.title.y = element_text(margin=margin(t=0, r=10.5, b=0, l=0), size=26),  #ditto
		axis.ticks = element_blank(),
		legend.position = c(.25,.9),
		legend.text = element_text(size = 19),
		legend.title=element_blank()
	)
