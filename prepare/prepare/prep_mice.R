source("./prepare/prep_env.R")
source("./prepare/functions.R")

# ---- load datasets ----
# mice global imputed episode data
df_ep_mice <- read.csv('./data/ep_uva_1B1A_mice.csv')
# mice global imputed safe window data
df_safe_window_mice <- read.csv('./data/safe48_allvars_uva_mice.csv')


# ---- keep positive patient and negative patient id ----
pos_id_lst <- unique(df_ep_mice$id[df_ep_mice$BSI==1])
neg_id_lst <- unique(df_ep_mice$id[df_ep_mice$BSI==0])
stopifnot(length(intersect(pos_id_lst, neg_id_lst)) == 0)


# ---- prepare dataframe subsets -----
df_ep_mice <- subset_df(df_ep_mice, neg_id_lst)
df_safe_window_mice <- subset_df(df_safe_window_mice, neg_id_lst)

