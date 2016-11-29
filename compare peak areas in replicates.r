library(gdata)
library(ggplot2)
library(stringr)
library(extrafont)

loadfonts()

setwd('C:\\Users\\Judson\\Documents\\plcgpaper\\All expt data dump--no protein name lookup')

sampleList = list.files()

tableOfR2 = c()

for(a in sampleList){
	for(b in sampleList){
		if(a==b){
			next
		}
		adat = read.xls(a)
		bdat = read.xls(b)

		aIDs = paste(adat$'assigned.sequence', adat$'MySQL.parsed.data..charge.state',sep='_')
		bIDs = paste(bdat$'assigned.sequence', bdat$'MySQL.parsed.data..charge.state',sep='_')

		adat$pepID = aIDs
		bdat$pepID = bIDs

		#parse filename & make experiment labels for plot
		aExptID = str_match(a, '(WT|KO)_IAP_([0-9]+min)_R([0-9]+)')
		bExptID = str_match(b, '(WT|KO)_IAP_([0-9]+min)_R([0-9]+)')
		#skip comparison if both replicates are same genotype
		if((aExptID[,2] != bExptID[,2]) | (aExptID[,3] != bExptID[,3])){
			next
		}

		aExptID = paste(aExptID[1,2]," ",aExptID[1,3],", replicate ",aExptID[1,4],sep="")
		bExptID = paste(bExptID[1,2]," ",bExptID[1,3],", replicate ",bExptID[1,4],sep="")
		aExptID = gsub("WT", "Jgamma1.WT", aExptID)
		aExptID = gsub("KO", "Jgamma1", aExptID)
		bExptID = gsub("WT", "Jgamma1.WT", bExptID)
		bExptID = gsub("KO", "Jgamma1", bExptID)

		#convert peak.RT and peak.area.final to numeric data type
		adat$peak.area.final <- as.numeric(as.character(adat$peak.area.final))
		bdat$peak.area.final <- as.numeric(as.character(bdat$peak.area.final))
		adat$peak.RT <- as.numeric(as.character(adat$peak.RT))
		bdat$peak.RT <- as.numeric(as.character(bdat$peak.RT))

		#Remove duplicate peptides (same label, peptide, and charge) in the data by selecting the observation with the highest peak area
		#(Must be disabled if there are no duplicates)
		newOrder <- order(adat$peak.area.final,decreasing=TRUE)
		adat <- adat[newOrder,]
		aIDs <- aIDs[newOrder]
		newOrder <- order(bdat$peak.area.final,decreasing=TRUE)
		bdat <- bdat[newOrder,]
		bIDs <- bIDs[newOrder]
		adat <- adat[-which(duplicated(aIDs)),]
		aIDs <- aIDs[-which(duplicated(aIDs))]
		bdat <- bdat[-which(duplicated(bIDs)),]
		bIDs <- bIDs[-which(duplicated(bIDs))]
		
		#Find match for each adat-detected peak in the bdat data
		dat <- data.frame(pep = adat$pepID[which(!is.na(match(aIDs,bIDs)))],
			  peakA = log( adat$peak.area.final[which(!is.na(match(aIDs,bIDs)))] ),
			  peakB = log( bdat$peak.area.final[na.omit(match(aIDs,bIDs))] ),
			  rtA = adat$peak.RT[which(!is.na(match(aIDs,bIDs)))],
			  rtB = bdat$peak.RT[na.omit(match(aIDs,bIDs))]
			)
	
		if(any(is.na(dat$peakA))){
			dropRows = which(is.na(dat$peakA))
			dat = dat[-dropRows,]
		}
		if(any(is.na(dat$peakB))){
			dropRows = which(is.na(dat$peakB))
			dat = dat[-dropRows,]
		}

		###Make the plot###
		#Different approach with ggplot. Reduce opacity of each point so that stacked datapoints appear darker
		r = range(c(dat$peakA, dat$peakB))
		r[1] = floor(r[1])
		r[2] = ceiling(r[2])

		peak.lm = lm(peakA ~ peakB, data=dat) 
		rsquared = summary(peak.lm)$r.squared

		# postscript(paste('C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\',a,' vs ',b,'.eps',sep=""), font = 'Arial Black')

		p = ggplot(dat,aes(x=peakA,y=peakB)) +
			geom_point(alpha = 0.5, shape = 16, col='navyblue', size=1) +
			ggtitle(bquote(R^2 == .(round(rsquared, 3)))) +
			xlab(bExptID) +
			ylab(aExptID) +
			scale_x_continuous(breaks=seq(r[1],r[2],by=4), limits=c(r[1],r[2]), labels=c(seq(r[1],r[2],by=4))) +
			scale_y_continuous(breaks=seq(r[1],r[2],by=4), limits=c(r[1],r[2]), labels=c(seq(r[1],r[2],by=4))) +
			stat_density2d(aes(fill = ..level..), geom = 'polygon', alpha = 0.015, bins = 150, colour = NA) +
			scale_fill_gradientn(colours=rev(rainbow(100, start=0, end=0.75))) +
			theme(
				title = element_text(color='black',size=32, family='Arial Black'),
				legend.title = element_blank(),
				legend.key = element_blank(),
				legend.text = element_text(color='black', size=24, family='Arial Black'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(), 
				panel.background = element_blank(),
				#axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
				axis.line.x = element_line(color='black', size = 0.8),
				axis.line.y = element_line(color='black', size = 0.8),
				axis.title.x = element_text(color='black', size = 24, family='Arial Black', vjust = 0.9),
				axis.title.y = element_text(color='black', size = 24, family='Arial Black', hjust = 0.9),
				axis.text = element_text(color = 'black', size = 24, family='Arial Black'),
				legend.position = 'none', #c(0.85,0.15),
				plot.margin=unit(c(t=15,r=17,b=17,l=17),'pt'),
				panel.border = element_rect(colour = "black", fill=NA, size=5)
			)
		ggsave(paste('C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\',a,' vs ',b,'.pdf',sep=''), plot = p)

		#Finally, add r2 to table
		tableOfR2 = rbind(tableOfR2, c(paste(aExptID,' vs. ',bExptID,sep=''), rsquared))
	}
}
#export R2 table
write.table(file='C:\\Users\\Judson\\Documents\\plcgpaper\\Supplemental material\\Rep comparison plots\\R2 for all replicate comparisons.txt',
  tableOfR2, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)