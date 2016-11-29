library(gdata)
library(ggplot2)
library(extrafont)

font_import()
setwd("C:/Users/Judson/Documents/plcgpaper/Line charts")
dat = read.xls('Jgamma1 all unique sites for line charts.xls')
#CD28 Y191 = 172
#CD28 Y206Y209 = 173
#CD28 Y218 = 174
#CD28 Y209 = 175
#Itk Y512 = 691
#Lck Y394 = 744
#Lck Y192 = 746
#Lck Y505 = 749
#Plcg2 Y753 = 990
#Plcg2 Y1217 = 991
#Plcg2 Y1245 = 992
#Plcg2 Y759 = 994
#TCRz Y142 = 1326
#TCRz Y153 = 1327
#TCRz Y111 = 1328
#TCRz Y123 = 1329
#TCRz S58Y72 = 1330
#TCRz S146Y153 = 1331
#TCRz Y142T149 = 1332
#TCRz Y142S146 = 1333
#TCRz Y142T147 = 1334
#TCRz Y64Y72 = 1335
#TCRz Y83 = 1336
#TCRz Y72 = 1337
#TCRz Y142T152 = 1338
#TCRz Y64 = 1339
#Zap70 Y493 = 1481
#Zap70 Y492 = 1484
#Zap70 Y292 = 1488
#Zap70 Y319 = 1493
#Zap70 Y315 = 1498

r = 691

peakDat = data.frame(timePoint=c(0,1,2,3,5,10), Jg1peakArea=rep(NA,6), Jg1peakSE=rep(NA,6), Jg1WTpeakArea=rep(NA,6), Jg1WTpeakSE=rep(NA,6))

#func to calculate std error of the mean
std <- function(x) sd(x)/sqrt(length(x))

#collate all peptide observations into a single row representing the protein they derive from
for(g in 2:3){
	qVals = c()
	for(tp in 1:6){
		peakAreas = c()
		for(rpt in 1:5){
			peakArea = dat[r,eval(paste('peakarea.manual.',g,'.rep',rpt,'.thresholded.timepoint',tp,sep=''))]
			if(is.na(peakArea) | peakArea < 1000){
				peakAreas = c(peakAreas, NA)
			} else {
				peakAreas = c(peakAreas, peakArea)
			}
		}
		if(g == 2){
			peakDat$Jg1WTpeakArea[tp] = mean(na.omit(peakAreas))
			peakDat$Jg1WTpeakSE[tp] = std(na.omit(peakAreas))
		} else if (g == 3){
			peakDat$Jg1peakArea[tp] = mean(na.omit(peakAreas))
			peakDat$Jg1peakSE[tp] = std(na.omit(peakAreas))
		}
		qVals = c(qVals, dat[r,eval(paste('qvalues.for.SILAC.timepoint',tp,sep=''))])
	}
}

#transform qvalues to * for representing significant data on the plot.
# assign * to the cell line with the higher peak area.
qVals = replace(qVals, which(qVals < .05), '*')
qVals = replace(qVals, which(qVals != '*'), '')
qVals = replace(qVals, which(is.na(qVals)), '')
qVals = c(qVals, qVals)

for(tp in 1:6){
	if(qVals[tp]=='*'){
		if(peakDat$Jg1WTpeakArea[tp] > peakDat$Jg1peakArea[tp]){
			qVals[tp+6] = ''
		} else if(peakDat$Jg1WTpeakArea[tp] < peakDat$Jg1peakArea[tp]){
			qVals[tp] = ''
		} else {
			print('WTF???')
		}
	}
}

###make peakDat2
peakDat2 = data.frame(timePoint=c(0,1,2,3,5,10,0,1,2,3,5,10),
  cellLine=c(rep('Jgamma1.WT',6),rep('Jgamma1',6)),
  peakArea=c(peakDat$Jg1WTpeakArea,peakDat$Jg1peakArea),
  peakStdErr = c(peakDat$Jg1WTpeakSE,peakDat$Jg1peakSE),
  sig = qVals)
####
r = range(peakDat2$peakArea)
r[2] = round(r[2], digits=4)*1.7

#postscript(fonts='Arial Black')

ggplot(peakDat2, aes(x=timePoint, y=peakArea)) +
  xlab('Time(minutes)') +
  ylab('Peak area') +
  geom_errorbar(aes(ymax=peakArea+peakStdErr, ymin=peakArea-peakStdErr, group=cellLine, color=cellLine), size=0.5, width=0.5) +
  geom_line(aes(y=peakArea, color=cellLine, group=cellLine), size=0.50) +
  scale_x_continuous(breaks=c(0,5,10)) +
  scale_y_continuous(breaks=c(0,r[2]), limits=c(0,r[2]), labels=c(0, format(r[2],scientific=TRUE, digits=1))) +
#  geom_point(data = peakDat2, aes(y=peakArea, shape=cellLine, group=cellLine), size=4.3) +
#  scale_shape_manual(values = c(2,16)) +
  geom_text(data = peakDat2, aes(x = timePoint, y = peakArea, label = sig), vjust = -1.4, size = 15) +
#  ggtitle(expression(paste('T cell receptor ', zeta, ' chain Y83'))) +
  ggtitle('Itk Y512') +
  theme(
    title = element_text(color='black',size=32),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(color='black', size=30),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    #axis.line = element_line(colour = 'black'), #doesn't work in all versions of ggplot2, so set lines individually instead
    axis.line.x = element_line(color='black', size = 0.8),
    axis.line.y = element_line(color='black', size = 0.8),
    axis.title.x = element_text(color='black', size = 30),
    axis.title.y = element_text(color='black', size = 30),
    axis.text = element_text(color = 'black', size = 29),
    legend.position = c(0.85,0.95),
    plot.margin=unit(c(t=1,r=75,b=12,l=5),'pt')
   )
