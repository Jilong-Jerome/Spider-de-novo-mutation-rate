library(tidyverse)

args <- commandArgs(T)
in_file <- args[1]
outname <- args[2]

region_df <- read_tsv(in_file)

region_df %>%
  mutate(end = end - 1)%>%
  mutate(len = end - start + 1)%>%
  mutate(lglen = log10(len))%>%
  filter(lglen > 2.5)%>%
  filter(poly_stat < 0.15)%>%
  select(region_n,start,end,poly_stat,len)%>%
  write_tsv(paste(outname,"tsv",sep="."))

png(paste(outname,"png",sep = "."))
ggplot()+
  geom_point(data = region_df%>%mutate(end = end - 1)%>%
  mutate(len = end - start + 1)%>%
  mutate(lglen = log10(len))%>%
  filter(lglen > 2.5)%>%
  filter(poly_stat < 0.15),mapping = aes(x=lglen,y=poly_stat),color = "black",alpha = 0.1)+
  geom_point(data = region_df%>%mutate(end = end - 1)%>%
  mutate(len = end - start + 1)%>%
  mutate(lglen = log10(len))%>%
  filter(lglen <= 2.5 | poly_stat >= 0.15),
  mapping = aes(x=lglen,y=poly_stat),color = "grey",alpha = 0.1)
dev.off()
