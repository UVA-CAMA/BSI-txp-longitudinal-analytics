source("./prepare/prep_env.R")
source("./prepare/functions.R")

# ---- load datasets ----
# median global imputed episode data
df_ep_median <- read.csv('./data/ep_uva_1B1A_median.csv')## 1B1A: 1 day before and 1 day after a culture
# median global imputed safe window data
df_safe_window_median <- read.csv('./data/safe48_allvars_uva_median.csv')

# ---- keep positive patient and negative patient id ----
pos_id_lst <- unique(df_ep_median$id[df_ep_median$BSI==1])
neg_id_lst <- unique(df_ep_median$id[df_ep_median$BSI==0])
stopifnot(length(intersect(pos_id_lst, neg_id_lst)) == 0)

# ---- subset dataframe ----
df_ep_median <- subset_df(df_ep_median, neg_id_lst)
df_safe_window_median <- subset_df(df_safe_window_median, neg_id_lst)
