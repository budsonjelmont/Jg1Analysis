library(gdata)
library(ggplot2)
library(stringr)
library(extrafont)

loadfonts()

setwd('C:\\Users\\Judson\\Documents\\plcgpaper\\All expt data dump--with labels')

#make empty data frame with same column structure to hold all collated runs
sampleList = list.files()
dat = read.xls(sampleList[1])
dat$pepID = NA
alldat = dat[FALSE,]

for(s in sampleList){
	dat = read.xls(s)
	IDs = paste(dat$'assigned.sequence', dat$'MySQL.parsed.data..charge.state',
	dat$'ExperimentGlobal..timecourse.position',sep='_')
	dat$pepID = IDs
	alldat = rbind(alldat, dat)
}

alldatWT = alldat[which(alldat$'ExperimentGlobal..timecourse.position' <= 6),]
alldatKO = alldat[which(alldat$'ExperimentGlobal..timecourse.position' >= 7),] 

r2table = c()

for(a in 1:5){
	for(b in 1:5){
		if(a==b){
			next
		}

		adatwt = alldatWT[which(alldatWT$'ExperimentGlobal..replicate.number'==a),]
		bdatwt = alldatWT[which(alldatWT$'ExperimentGlobal..replicate.number'==b),]
		adatko = alldatKO[which(alldatKO$'ExperimentGlobal..replicate.number'==a),]
		bdatko = alldatKO[which(alldatKO$'ExperimentGlobal..replicate.number'==b),]

		#parse filename & make experiment labels for plot
		awtExptID = paste("Jgamma1.WT rep",a,sep="")
		bwtExptID = paste("Jgamma1.WT rep",b,sep="")
		akoExptID = paste("Jgamma1 rep",a,sep="")
		bkoExptID = paste("Jgamma1 rep",b,sep="")

		#convert peak.RT and peak.area.final to numeric data type
		adatwt$peak.area.final <- as.numeric(as.character(adatwt$peak.area.final))
		bdatwt$peak.area.final <- as.numeric(as.character(bdatwt$peak.area.final))
		adatwt$peak.RT <- as.numeric(as.character(adatwt$peak.RT))
		bdatwt$peak.RT <- as.numeric(as.character(bdatwt$peak.RT))
		adatko$peak.area.final <- as.numeric(as.character(adatko$peak.area.final))
		bdatko$peak.area.final <- as.numeric(as.character(bdatko$peak.area.final))
		adatko$peak.RT <- as.numeric(as.character(adatko$peak.RT))
		bdatko$peak.RT <- as.numeric(as.character(bdatko$peak.RT))

		#Remove duplicate peptides (same label, peptide, and charge) in the data by selecting the observation with the highest peak area
		#(Must be disabled if there are no duplicates)
		awtIDs = adatwt$pepID
		bwtIDs = bdatwt$pepID
		akoIDs = adatko$pepID
		bkoIDs = bdatko$pepID
		
		newOrder <- order(adatwt$peak.area.final,decreasing=TRUE)
		adatwt <- adatwt[newOrder,]
		awtIDs <- awtIDs[newOrder]
		newOrder <- order(bdatwt$peak.area.final,decreasing=TRUE)
		bdatwt <- bdatwt[newOrder,]
		bwtIDs <- bwtIDs[newOrder]
		newOrder <- order(adatko$peak.area.final,decreasing=TRUE)
		adatko <- adatko[newOrder,]
		akoIDs <- akoIDs[newOrder]
		newOrder <- order(bdatko$peak.area.final,decreasing=TRUE)
		bdatko <- bdatko[newOrder,]
		bkoIDs <- bkoIDs[newOrder]
		
		if(any(duplicated(awtIDs))){
			adatwt <- adatwt[-which(duplicated(awtIDs)),]
			awtIDs <- awtIDs[-which(duplicated(awtIDs))]
			}
		if(any(duplicated(bwtIDs))){
			bdat <- bdatwt[-which(duplicated(bwtIDs)),]
			bwtIDs <- bwtIDs[-which(duplicated(bwtIDs))]
		}
		if(any(duplicated(akoIDs))){
			adatko <- adatko[-which(duplicated(akoIDs)),]
			akoIDs <- akoIDs[-which(duplicated(akoIDs))]
			}
		if(any(duplicated(bkoIDs))){
			bdat <- bdatko[-which(duplicated(bkoIDs)),]
			bkoIDs <- bkoIDs[-which(duplicated(bkoIDs))]
		}
		#Find match for each adat-detected peak in the bdatwt data
		datwt <- data.frame(pep = adatwt$pepID[which(!is.na(match(awtIDs,bwtIDs)))],
			  peakA = log( adatwt$peak.area.final[which(!is.na(match(awtIDs,bwtIDs)))] ),
			  peakB = log( bdatwt$peak.area.final[na.omit(match(awtIDs,bwtIDs))] ),
			  rtA = adatwt$peak.RT[which(!is.na(match(awtIDs,bwtIDs)))],
			  rtB = bdatwt$peak.RT[na.omit(match(awtIDs,bwtIDs))]
			)
		datko <- data.frame(pep = adatko$pepID[which(!is.na(match(akoIDs,bkoIDs)))],
			  peakA = log( adatko$peak.area.final[which(!is.na(match(akoIDs,bkoIDs)))] ),
			  peakB = log( bdatko$peak.area.final[na.omit(match(akoIDs,bkoIDs))] ),
			  rtA = adatko$peak.RT[which(!is.na(match(akoIDs,bkoIDs)))],
			  rtB = bdatko$peak.RT[na.omit(match(akoIDs,bkoIDs))]
			)

		if(any(is.na(datwt$peakA))){
			dropRows = which(is.na(datwt$peakA))
			datwt = datwt[-dropRows,]
		}
		if(any(is.na(datwt$peakB))){
			dropRows = which(is.na(datwt$peakB))
			datwt = datwt[-dropRows,]
		}
		if(any(is.na(datko$peakA))){
			dropRows = which(is.na(datko$peakA))
			datko = datko[-dropRows,]
		}
		if(any(is.na(datko$peakB))){
			dropRows = which(is.na(datko$peakB))
			datko = datko[-dropRows,]
		}

		###Make the plot###
		###WT first
		
		#override dynamic axis ranges
		# r = range(c(datwt$peakA, datwt$peakB))
		# r[1] = floor(r[1])
		# r[2] = ceiling(r[2])
		r = c(9,26)

		peak.lm = lm(peakA ~ peakB, data=datwt) 
		rsquared = summary(peak.lm)$r.squared

		# postscript(paste('C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\',a,' vs ',b,'.eps',sep=""), font = 'Arial')

		p = ggplot(datwt,aes(x=peakA,y=peakB)) +
			geom_point(alpha = 0.5, shape = 16, col='navyblue', size=1) +
			ggtitle(bquote(R^2 == .(round(rsquared, 3)))) +
			xlab(bwtExptID) +
			ylab(awtExptID) +
			scale_x_continuous(breaks=seq(r[1],r[2],by=4), limits=c(r[1],r[2]), labels=c(seq(r[1],r[2],by=4))) +
			scale_y_continuous(breaks=seq(r[1],r[2],by=4), limits=c(r[1],r[2]), labels=c(seq(r[1],r[2],by=4))) +
			stat_density2d(aes(fill = ..level..), geom = 'polygon', alpha = 0.015, bins = 150, colour = NA) +
			scale_fill_gradientn(colours=rev(rainbow(100, start=0, end=0.75))) +
			theme(
				title = element_text(color='black',size=32, family='Arial'),
				legend.title = element_blank(),
				legend.key = element_blank(),
				legend.text = element_text(color='black', size=24, family='Arial'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), 
				panel.background = element_blank(),
				#axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
				axis.line.x = element_line(color='black', size = 0.8),
				axis.line.y = element_line(color='black', size = 0.8),
				axis.title.x = element_text(color='black', size = 26, family='Arial', vjust = 0.9, margin = unit(c(t=4, r=0, b=0, l=0), 'mm')),
				axis.title.y = element_text(color='black', size = 26, family='Arial', hjust = 0.5, vjust=0.99, margin = unit(c(t=3, r=4, b=0, l=0), 'mm')),
				axis.text = element_text(color = 'black', size = 25, family='Arial'),
				legend.position = 'none', #c(0.85,0.15),
				plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
				panel.border = element_rect(colour = "black", fill=NA, size=1.5)
			)
		ggsave(paste('C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\Jgamma1.WT ',a,' vs ',b,'.pdf',sep=''), plot = p)

		#Finally, add r2 to table
		r2table = rbind(r2table, c(paste(awtExptID,' vs. ',bwtExptID,sep=''), rsquared))
		
		
		###KO
		# r = range(c(datko$peakA, datwt$peakB))
		# r[1] = floor(r[1])
		# r[2] = ceiling(r[2])

		peak.lm = lm(peakA ~ peakB, data=datko) 
		rsquared = summary(peak.lm)$r.squared

		# postscript(paste('C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\',a,' vs ',b,'.eps',sep=""), font = 'Arial')

		p = ggplot(datko,aes(x=peakA,y=peakB)) +
			geom_point(alpha = 0.5, shape = 16, col='navyblue', size=1) +
			ggtitle(bquote(R^2 == .(round(rsquared, 3)))) +
			xlab(bkoExptID) +
			ylab(akoExptID) +
			scale_x_continuous(breaks=seq(r[1],r[2],by=4), limits=c(r[1],r[2]), labels=c(seq(r[1],r[2],by=4))) +
			scale_y_continuous(breaks=seq(r[1],r[2],by=4), limits=c(r[1],r[2]), labels=c(seq(r[1],r[2],by=4))) +
			stat_density2d(aes(fill = ..level..), geom = 'polygon', alpha = 0.015, bins = 150, colour = NA) +
			scale_fill_gradientn(colours=rev(rainbow(100, start=0, end=0.75))) +
			theme(
				title = element_text(color='black',size=32, family='Arial'),
				legend.title = element_blank(),
				legend.key = element_blank(),
				legend.text = element_text(color='black', size=24, family='Arial'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), 
				panel.background = element_blank(),
				#axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
				axis.line.x = element_line(color='black', size = 0.8),
				axis.line.y = element_line(color='black', size = 0.8),
				axis.title.x = element_text(color='black', size = 26, family='Arial', vjust = 0.9, margin = unit(c(t=4, r=0, b=0, l=0), 'mm')),
				axis.title.y = element_text(color='black', size = 26, family='Arial', hjust = 0.5, vjust=0.99, margin = unit(c(t=3, r=4, b=0, l=0), 'mm')),
				axis.text = element_text(color = 'black', size = 25, family='Arial'),
				legend.position = 'none', #c(0.85,0.15),
				plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
				panel.border = element_rect(colour = "black", fill=NA, size=1.5)
			)
		ggsave(paste('C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\Jgamma1 ',a,' vs ',b,'.pdf',sep=''), plot = p)

		#Finally, add r2 to table
		r2table = rbind(r2table, c(paste(akoExptID,' vs. ',bkoExptID,sep=''), rsquared))
	}
}
#export R2 table
write.table(file='C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\R2 for all replicate comparisons.txt',
  r2table, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)