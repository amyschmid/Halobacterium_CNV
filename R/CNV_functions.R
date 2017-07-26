############################################################
##
##  Functions for analysis of Halobacterium salinarum
##  CNV/microarray data
##  Keely Dulmage
##  Schmid Lab, Duke University
##
############################################################

#threshold both the amplification/deletion threshold: "threshold", and the size of the fragment "f.size" 
copymap4=function(ref_array, seg_array,threshold,f.size){
  for(i in 1:length(ref_array[,1])){
    row_ident=subset(seg_array,(seg_array$chrom==ref_array$Chr[i]) & (seg_array$loc.start<=ref_array$Start[i]) & (seg_array$loc.end>= ref_array$Start[i]) & ((seg_array$loc.end-seg_array$loc.start)>=f.size))
    numb_up=length(row_ident$seg.mean[row_ident$seg.mean>=0])
    numb_down=length(row_ident$seg.mean[row_ident$seg.mean<=0])
    thresh_up=length(row_ident$seg.mean[row_ident$seg.mean>=threshold])
    thresh_down=length(row_ident$seg.mean[row_ident$seg.mean<=-(threshold)])
    av_ratio=mean(row_ident$seg.mean)
    ref_array[i,4]=numb_up
    ref_array[i,5]=numb_down
    ref_array[i,6]=thresh_up
    ref_array[i,7]=thresh_down
    ref_array[i,8]=av_ratio
  }
  ref_array
} 

#Split checks for completeness (no NAs), orders the table by locus, and splits it into three tables by chromosome
#Each chromosome is accessed like so: "split.map[[3]]" would give the third chromosome

split2=function(copymap){
  copymap2=na.omit(copymap) 
  copy.s=copymap2[order(copymap2$Start),]
  copy=copy.s[order(copy.s$Chr),]
  Chr=subset(copy,copy$Chr=='Chr')
  pNRC100=subset(copy,copy$Chr=='pNRC100')
  pNRC200=subset(copy,copy$Chr=='pNRC200')
  array_list=list(Chr,pNRC100,pNRC200)
  array_list } 


#plot chromosome for chip
plotchr=function(copymap,narrays,ymax){
  all_thresh=copymap$Thresh_up + copymap$Threshdown
  plot(copymap$Start,((all_thresh/narrays)*100),xlab="Genome Position",ylab='Frequency (%)',main='', ylim = c(0,ymax),pch=16, cex.axis=1.5,cex.lab=1.5,type='l', lwd=1.5)
}

add_copyplot=function(copymap,narrays,line_color){
  all_thresh=copymap$Thresh_up + copymap$Threshdown
  points(copymap$Start,((all_thresh/narrays)*100),col=line_color,type='l',lty=2, lwd=1.5)}


plotplas.GE2=function(copymap,narrays,yrange){
  plot(copymap$Start,((copymap$Thresh_up/narrays)*100),xlab="Genome Position",ylab='Frequency (%)',main='', ylim = yrange, pch=16, cex.axis=1.5,cex.lab=1.5,col='orange',type='l',lwd=1.5)
  points(copymap$Start,((copymap$Threshdown/narrays)*100), type='l', lwd=1.5, col='blue')
  abline(h=0,col='gray')
  #abline(v=1806766,lty=1, col='red',lty=2)
  legend('topright',legend=c('Upregulation','Downregulation'),col=c('orange','blue'),bg='white', lty=c(1,1),lwd=2)
}


copyplot.chr=function(segmented_object,ylim_range){
  datatest2=data.frame(segmented_object$data[[1]],segmented_object$data[[2]],segmented_object$data[[3]])
  colnames(datatest2)=c('Chr','Start','Intensity')
  datatest.chr=datatest2[datatest2$Chr=='Chr',]
  #sorts segments by chromosome
  allseg=segmented_object[[2]]
  allseg.chr=allseg[allseg$chrom=='Chr',]
  
  #plot chromosome and segments
  plot(datatest.chr$Start,datatest.chr$Intensity,pch=16,ylim=ylim_range, cex=.4, xlab="Genome Position (bp)",main='Chromosome',cex.lab=1.5,cex.axis=1.5,ylab='',cex.main=2)
  abline(h=0,col='lightgray',lwd=3)
  abline(h=c(-0.5, 0.5), col = "purple", lty = 2)
  legend('topright',legend=c('Segment','Threshold'),col=c('green','purple'),bg='white', lty=c(1,2),lwd=2)
  for(i in 1:length(allseg.chr[,1])){
    segments(allseg.chr$loc.start[i],allseg.chr$seg.mean[i], allseg.chr$loc.end[i],allseg.chr$seg.mean[i],col='green',lwd=4)
  }		
}


#Plots amplification on pos Y axis, deletion on neg Y axis
plotchr_chip7=function(copymap,narrays,xrange, ish_table){
  plot(copymap$Start,((copymap$Thresh_up/narrays)*100),xlab="Genome Position",ylab='% arrays',main='Chromosome, frequency of copy number events', ylim = c(-25,75),xlim=xrange,pch=16, cex.axis=1.5,cex.lab=1.5,col='orange',type='l',lwd=1.5)
  points(copymap$Start,((-copymap$Threshdown/narrays)*100), type='l', lwd=1.5, col='blue')
  #abline(v=920669,lty=1, col='red')
  #abline(v=1806766,lty=1, col='red')
  #abline(v=ish_table$Start,lty=2, col=rgb(0.1,0.1,.8,.3),lwd=2)
}


plotchr_chip2=function(copymap,narrays,ymax, ish_table){
  all_thresh=copymap$Thresh_up + copymap$Threshdown
  plot(copymap$Start,((all_thresh/narrays)*100),xlab="Genome Position",ylab='% arrays',main='Chromosome, frequency of copy number events', ylim = c(0,ymax),pch=16, cex.axis=1.5,cex.lab=1.5,col='black',type='l',lwd=1.5)
  #abline(v=920669,lty=1, col='red')
  #abline(v=1806766,lty=1, col='red')
  abline(v=ish_table$Start,lty=2, col=rgb(0.1,0.1,.8,.3),lwd=1.5)
  points(copymap$Start,((all_thresh/narrays)*100),col='black',type='l',lwd=1.5)
}


