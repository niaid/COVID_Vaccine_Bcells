
new_tp <- function(tp_vec){
  
  timepoint_order_orig <- c(
    "D1_pre",
    "D1_1week",
    "D1_10days",
    "D1_2week",
    "D2_pre",
    "D2_5days",
    "D2_1week",
    "D2_10days",
    "D2_2weeks",
    "D2_4weeks"
  )
  
  timepoint_order_new <- c(
    "v1D0",
    "v1D7",
    "v1D10",
    "v1D14",
    "v2D0",
    "v2D5",
    "v2D7",
    "v2D10",
    "v2D14",
    "v2D28"
  )
  
  out <- timepoint_order_new[match(tp_vec, timepoint_order_orig)]
  
}
