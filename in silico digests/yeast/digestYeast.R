### modified from R-bloggers.com
### load necessary libraries
library(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer3)  ## use genome of interest
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)

# Identify XhoI recognition sites for each chromosomal entry
# Generate dataframe with length of XhoI digested fragments 
mdf = data.frame();
for (i in seq_along(Scerevisiae)){
  print(paste("Processing ", seqnames(Scerevisiae)[i], sep = ""))
  m <- matchPattern("CTCGAG", Scerevisiae[[i]]) ## insert your enzyme's restriction site!!
  
  # try to do a double digest using n
  n <- matchPattern("GAATTC", m)
  starts <- start(gaps(n))
  ends <- end(gaps(n))
  temp_df <- data.frame(start = starts-4, end = ends, chr = seqnames(Scerevisiae)[i])
  temp_df$start <- replace(temp_df$start, temp_df$start == -3, 0)
  temp_df <- temp_df[c("chr", "start", "end")]
  mdf <- rbind(mdf, temp_df)
}

# Extract digested fragment length
# Plot fragments within desired size range

mdf$width = mdf$end - mdf$start
ml <- mdf[mdf$width>0 & mdf$width<601,]
counts <- ddply(ml, .(width), nrow)

# plot the frequency of the fragment lengths (y-axis is logarithmic)

plot <- ggplot(counts, aes(x=width, y=V1)) + geom_line()
plot + scale_y_continuous(trans = log2_trans())

# write the counts to a .csv file
write.table(counts,"yeasttest1.csv",sep="\t",row.names=FALSE)
