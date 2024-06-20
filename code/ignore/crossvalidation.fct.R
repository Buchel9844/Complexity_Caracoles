k_split <- function(dat,k, seed = 1){
  set.seed(seed)
  grp <- rep(1:k, length.out = nrow(dat))
  dat |>
    mutate(grp = sample(grp, nrow(dat), replace = F)) |>
    group_split(grp)|>
    map(\(d) select(d, -grp))
}

SpData <- SpData[which(SpData$year == as.numeric(year)),] %>%
  dplyr::filter(focal == FocalPrefix)

list.StData <- k_split(SpData,5)


  
  SpData <- list.StData[   c(1:5)[!1:5 %in% split.data]] %>%
    bind_rows()
  
  test <- list.StData[[split.data]]