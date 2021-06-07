source("./prepare/prep_env.R")
source("./prepare/functions.R")

# ---- load datasets ----
# no global imputed episode data (_org = no global imputation yes (12+24) ffill)
df_ep_raw <- read.csv('./data/ep_uva_1B2A_raw.csv') # 
# no global imputed safe window data (_org = no global imputation yes (12+24) ffill)
df_safe_window_raw <- read.csv('./data/safe48_allvars_uva_raw.csv')


# ---- keep positive patient and negative patient id ----
pos_id_lst <- unique(df_ep_raw$id[df_ep_raw$BSI==1])
neg_id_lst <- unique(df_ep_raw$id[df_ep_raw$BSI==0])
stopifnot(length(intersect(pos_id_lst, neg_id_lst)) == 0)

# ---- prepare dataframe subset ----
df_ep_raw <- subset_df(df_ep_raw, neg_id_lst)
df_safe_window_raw <- subset_df(df_safe_window_raw, neg_id_lst)
