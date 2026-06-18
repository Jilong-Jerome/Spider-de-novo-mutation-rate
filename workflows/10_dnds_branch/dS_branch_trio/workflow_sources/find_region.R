library(tidyverse)

args <- commandArgs(T)
in_file <- args[1]
outname <- args[2]
print(in_file)
print(outname)
mask_dt <- read_tsv(in_file)

mask_df <- mask_dt %>%
  mutate(last = lag(remove))%>%
  mutate(bound = ifelse(last == TRUE & remove == FALSE,"start",
                        ifelse(last == FALSE & remove == TRUE, "end", "non")))
if (mask_dt$remove[1] == FALSE){
  mask_df$bound[1] = "start"
}
if (mask_dt$remove[1] == TRUE){
  mask_df$bound[1] = "non"
}
if (mask_dt$remove[length(mask_dt$remove)] == FALSE){
  temp <- mask_df[nrow(mask_df),]
  mask_df[nrow(mask_df)+1,] <- temp
  mask_df$bound[nrow(mask_df)] = "end"
  # Use exclusive end (site+1) to match mid-region convention,
  # so filter_region.R's end-1 produces the correct inclusive boundary
  mask_df$site[nrow(mask_df)] = mask_df$site[nrow(mask_df)] + 1
}

bound_dt <- mask_df%>%filter(bound != "non")
re_n = 1
re_start_vec <- c()
re_n_vec <- c()
re_end_vec <- c()
for (i in 1:length(bound_dt$bound)){
  if (bound_dt$bound[i] == "start"){
    re_n_vec <- c(re_n_vec,re_n)
    re_start_vec <- c(re_start_vec,bound_dt$site[i])
  }
  if (bound_dt$bound[i] == "end"){
    re_end_vec <- c(re_end_vec,bound_dt$site[i])
    re_n <- re_n + 1
  }
}

region_df <- data.frame(region_n = re_n_vec,
                        start = re_start_vec,
                        end = re_end_vec)
poly_vec <- c()
for (i in 1:length(region_df$region_n)){
  start <-region_df$start[i]
  end <- region_df$end[i]
  poly_num <- mask_df%>%
    filter(site >= start)%>%
    filter(site <= end)%>%
    filter(poly == TRUE)%>%
    NROW()
  poly_stat <- poly_num/(end - start + 1)
  poly_vec <- c(poly_vec,poly_stat)
}
region_df$poly_stat <- poly_vec
write_tsv(region_df,outname)
