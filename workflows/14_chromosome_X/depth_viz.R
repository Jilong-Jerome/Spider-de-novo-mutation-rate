library(tidyverse)
args <- commandArgs(T)
# Create a custom theme
my_theme <- function() {
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),    # Remove major grid lines
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Add x and y axis lines
    legend.position = "right"
  )
}
theme_set(my_theme())

depth_file <- args[1]
depth_df <- read_tsv(depth_file,col_names = c("chrom","pos","depth"))

p <- depth_df %>%
	group_by(chrom)%>%
	summarise(median = median(depth),
		  q975 = quantile(depth,0.975),
		  q025 = quantile(depth,0.025),
		  mean = mean(depth))%>%
	ggplot()+
	geom_segment(aes(x= chrom,xend = chrom,y=q025,yend = q975))+
	geom_point(aes(x=chrom,y=median,color = "median"))+
	geom_point(aes(x=chrom,y=mean,color = "mean"))+
	xlab("chromosome")+
	ylab("depth")+
	scale_colour_manual(name = "Statistics",
                      values = c("median" = "black", "mean" = "red"),
                      labels = c("Median", "Mean"))
outname <- sys.argv[2]
ggsave(paste(outname,".png",sep=""),p,height = 4, width = 12)
