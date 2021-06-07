source("./prepare/prep_env.R")
source("./prepare/functions.R")

# ---- load datasets ----
# no global imputed episode data (_org = no global imputation yes (12+24) ffill)
df_ep_original <- read.csv('./data/ep_uva_1B2A_org.csv') # 
# no global imputed safe window data (_org = no global imputation yes (12+24) ffill)
df_safe_window_original <- read.csv('./data/safe48_allvars_uva_org.csv')

# ---- keep positive patient and negative patient id ----
pos_id_lst <- unique(df_ep_original$id[df_ep_original$BSI==1])
neg_id_lst <- unique(df_ep_original$id[df_ep_original$BSI==0])
stopifnot(length(intersect(pos_id_lst, neg_id_lst)) == 0)

# ---- prepare dataframe subset ----
df_ep_original <- subset_df(df_ep_original, neg_id_lst)
df_safe_window_original <- subset_df(df_safe_window_original, neg_id_lst)

# fulfill age in df_ep_original
for (ID in unique(df_ep_original$id)){
  df_ep_original$age[df_ep_original$id == ID] <- round(mean(df_ep_original$age[df_ep_original$id == ID], na.rm = TRUE))
}


