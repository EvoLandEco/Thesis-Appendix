tree_list <- list()

# Grab all chronograms
for (i in 1:length(study_list)) {
  tree_list[[length(tree_list) + 1]] <- rotl::get_study(study_id = study_list[i], object_format = "phylo")
}
tree_list_bak <- tree_list

# Data cleaning
cleaned_list <- list()
for (i in 1:length(tree_list)) {
  if (length(tree_list[[i]]) == 3 && inherits(tree_list[[i]], "phylo")) {
    next
  } else if (length(tree_list[[i]]) == 4 && inherits(tree_list[[i]], "phylo")) {
    cleaned_list[[length(cleaned_list) + 1]] <- tree_list[[i]]
  } else {
    for (j in 1:length(tree_list[[i]])) {
      if (length(tree_list[[i]][[j]]) == 3 && inherits(tree_list[[i]][[j]], "phylo")) {
        next
      } else if (length(tree_list[[i]][[j]]) == 4 && inherits(tree_list[[i]][[j]], "phylo")) {
        cleaned_list[[length(cleaned_list) + 1]] <- tree_list[[i]][[j]]
      } else {
        print(paste("Unexpected list structure in tree ", i, " and sublist ", j))
        stop("Unexpected list structure")
      }
    }
  }
}

# Check data quality
all_phylo <- all(sapply(cleaned_list, function(x) inherits(x, "phylo")))
all_time_calibrated <- all(sapply(cleaned_list, function(x) identical(names(x), c("edge", "tip.label", "Nnode", "edge.length"))))
print(all_phylo && all_time_calibrated)

# Remove very large trees
filtered_list <- Filter(function(x) x$Nnode < 1001, cleaned_list)

all_binary <- all(sapply(filtered_list, function(x) is.binary(x)))

filtered_list2 <- Filter(function(x) is.binary(x), filtered_list)

all_binary <- all(sapply(filtered_list2, function(x) is.binary(x)))

filtered_list3 <- Filter(function(x) is.ultrametric(x, tol = 0.00001), ott_trees)

filtered_list3 <- Filter(function(x) x$Nnode > 10, filtered_list3)

filtered_list3 <- sapply(filtered_list3, function(x) phytools::force.ultrametric(x), simplify = F)

filtered_list3 <- Filter(function(x) is.ultrametric(x), filtered_list3)

all_ultrametric <- all(sapply(filtered_list3, function(x) is.ultrametric(x)))

hist(sapply(filtered_list3, function(x) x$Nnode * 2 + 1), breaks = 20)

saveRDS(filtered_list3, "c:/Source/eveGNN/emp/OTT_trees.rds")

ott_mle_result <- load_empirical_ott_mle_result("D:/Habrok/Data/MLE")

ott_mle_result_ddd <- ott_mle_result$DDD

ott_gnn_result_ddd <- readRDS("D:/Habrok/Data/STBO/empirical_gnn_1_lstm_result.rds")

ott_gnn_result_ddd <- ott_gnn_result_ddd %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000) %>%
  dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                         num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                         num_nodes >= 500 ~ "Large"))

joint_ott_result <- dplyr::left_join(ott_mle_result_ddd, ott_gnn_result_ddd, by = "index")

emp_p1 <- ggplot(joint_ott_result) + geom_point(aes((lambda - mu), (pred_lambda_after - pred_mu_after), color = small)) +
  xlim(0,3) + ylim(0,3)  +
  geom_abline(slope = 1, color = "red", lty="dashed") +
  labs(x="MLE Prediction", y = "NN Prediction", title = "Net diversification rate") +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
     legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))

emp_p2 <- ggplot(joint_ott_result) + geom_point(aes(cap, pred_cap_after, color = small)) +
  xlim(0,900) + ylim(0,900) +
  geom_abline(slope = 1, color = "red", lty="dashed") +
  labs(x="MLE Prediction", y = "NN Prediction", title = "Carrying capacity") +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))

emp_p3 <- ggplot(joint_ott_result) + geom_point(aes(pred_cap_before, pred_cap_after, color = small)) +
  xlim(0,900) + ylim(0,900) +
  geom_abline(slope = 1, color = "red", lty="dashed") +
  labs(x="GNN Prediction", y = "BoostBT Prediction", title = "Carrying capacity") +
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))

ggthemr::ggthemr("flat")

emp_p1 + emp_p2 + plot_layout(nrow = 1, guides = "collect")

PBD::pbd_mean_durspecs(eris = c(0.2, 0.2), scrs = 1, siris = c(0.8, 0.6))

params <- yaml::read_yaml("../Config/ddd_sim.yaml")

if (!dir.exists(name)) {
  dir.create(name)
}

setwd(name)

dists <- params$dists
cap_range <- params$cap_range
max_mu <- params$max_mu
within_ranges <- params$within_ranges
nrep <- 50000
age <- params$age
ddmodel <- params$ddmodel
nworkers_sim <- params$nworkers_sim

future::plan("multisession", workers = 8)

ddd_free_tes_list <- future.apply::future_replicate(nrep, eveGNN::randomized_ddd_fixed_age(dists,
                                                                                           cap_range = cap_range,
                                                                                           max_mu = max_mu,
                                                                                           age = age,
                                                                                           model = ddmodel), simplify = FALSE)

ddd_list_all <- purrr::transpose(ddd_free_tes_list)

ddd_list_all <- readRDS("D:\\Habrok\\Data\\MLE\\DDD_MLE_TES\\MLE_DATA\\ddd_mle.rds")
future::plan("multisession", workers = 4)
ddd_summary_stats <- furrr::future_map(ddd_list_all$tes, eveGNN:::tree_to_stats)

# Flatten the lists of pars to a data frame
pars_cca <- ddd_list_all$pars %>% purrr::transpose()
names(pars_cca) <- c("lambda", "mu", "K")
pars_cca <- data.frame(lapply(pars_cca, function(sublist) {
  unlist(sublist, use.names = FALSE)
}))

# Flatten the lists of stats to a data frame
stats_cca <- ddd_summary_stats %>% dplyr::bind_rows()

write.csv2(pars_cca, "pars_cca.csv")
write.csv2(stats_cca, "stats_cca.csv")

ggthemr::ggthemr("flat", layout = "minimal")

df_beta <- data.frame(size = 2000* rbeta(100000,0.1,2))

df_unif <- data.frame(size = rep(1:2000, each = 50))

ggplot(df_beta) + geom_histogram(aes(size),bins=10)
ggplot(df_unif) + geom_histogram(aes(size),bins=20)


test_poly_gnn <- load_final_difference_by_layer("D:\\Habrok\\Data\\POLY", task_type = "DDD_POLY_TES", model_type = "diffpool", depth = 2)

test_poly_gnn2 <- test_poly_gnn %>% mutate_at(vars(contains("cap")),  ~ round(. * 1000, 0)) %>%
  group_by(lambda, .add = T) %>%
  group_by(mu, .add = T) %>%
  group_by(cap, .add = T) %>%
  summarise(mae_lambda = mean(abs(lambda-lambda_pred)),
            mae_mu = mean(abs(mu-mu_pred)),
            mae_cap = mean(abs(cap-cap_pred)),
            mean_node = mean(num_nodes)) %>%
  mutate(cap = if_else(cap == 700, true = 800, false = cap))

saveRDS(test_poly_gnn2, "D:\\Habrok\\Data\\POLY\\robustness_gnn.rds")

test_poly_dnn_lstm <- readRDS("D:\\Habrok\\Data\\POLY\\DDD_POLY_TES\\DDD_POLY_TES_dnn_lstm_1.rds")

test_poly_dnn_lstm2 <- test_poly_dnn_lstm %>% mutate_all(as.numeric) %>%
  mutate_at(vars(contains("cap")),  ~ round(. * 1000, 0)) %>%
  group_by(lambda, .add = T) %>%
  group_by(mu, .add = T) %>%
  group_by(cap, .add = T) %>%
  summarise(mae_lambda_dnn = mean(abs(lambda-lambda_pred_dnn)),
            mae_mu_dnn = mean(abs(mu-mu_pred_dnn)),
            mae_cap_dnn = mean(abs(cap-cap_pred_dnn)),
            mae_lambda_lstm = mean(abs(lambda-lambda_pred_lstm)),
            mae_mu_lstm = mean(abs(mu-mu_pred_lstm)),
            mae_cap_lstm = mean(abs(cap-cap_pred_lstm))) %>%
  mutate(cap = if_else(cap == 700, true = 800, false = cap))

test_poly_lstm <- readRDS("D:\\Habrok\\Data\\POLY\\DDD_POLY_TES\\DDD_POLY_TES_gnn_2_lstm_compensation.rds")

test_poly_lstm2 <- test_poly_lstm %>% mutate_at(vars(contains("cap")),  ~ round(. * 1000, 0)) %>%
  group_by(lambda, .add = T) %>%
  group_by(mu, .add = T) %>%
  group_by(cap, .add = T) %>%
  summarise(mae_lambda = mean(abs(lambda-pred_lambda_after)),
            mae_mu = mean(abs(mu-pred_mu_after)),
            mae_cap = mean(abs(cap-pred_cap_after))) %>%
  mutate(cap = if_else(cap == 700, true = 800, false = cap))

saveRDS(test_poly_lstm2, "C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\robustness\\robustness_Boost_BT.rds")

test_poly_lstm_dnn <- readRDS("D:\\Habrok\\Data\\POLY\\DDD_POLY_TES\\DDD_POLY_TES_gnn_2_lstm_compensation_after_dnn.rds")

test_poly_lstm_dnn2 <- test_poly_lstm_dnn %>% mutate_at(vars(contains("cap")),  ~ round(. * 1000, 0)) %>%
  group_by(lambda, .add = T) %>%
  group_by(mu, .add = T) %>%
  group_by(cap, .add = T) %>%
  summarise(mae_lambda = mean(abs(lambda-pred_lambda_after_lstm_after_dnn)),
            mae_mu = mean(abs(mu-pred_mu_after_lstm_after_dnn)),
            mae_cap = mean(abs(cap-pred_cap_after_lstm_after_dnn))) %>%
  mutate(cap = if_else(cap == 700, true = 800, false = cap))

saveRDS(test_poly_lstm_dnn2, "C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\robustness\\robustness_Boost_SS+BT.rds")

test_poly_mle_opt <- load_separated_mle_result(path = "D:\\Habrok\\Data\\POLY", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

test_poly_mle_opt2 <- test_poly_mle_opt %>% dplyr::filter(lambda_r_diff < 10000, lambda_r_diff>-10000, mu_r_diff<10000,mu_r_diff>-10000, cap_r_diff<10000, cap_r_diff>-10000)%>%
  group_by(lambda, .add = T) %>%
  group_by(mu, .add = T) %>%
  group_by(cap, .add = T) %>%
  summarise(mae_lambda = mean(abs(lambda_a_diff)),
            mae_mu = mean(abs(mu_a_diff)),
            mae_cap = mean(abs(cap_a_diff)),
            mean_node = mean(num_nodes)) %>%
  mutate(cap = if_else(cap == 700, true = 800, false = cap))

saveRDS(test_poly_mle_opt2, "C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\robustness\\robustness_MLE_best.rds")

test_poly_mle <- load_separated_mle_result(path = "D:\\Habrok\\Data\\POLY", task_type = "DDD", model_type = "diffpool", no_init = TRUE)

test_poly_mle2 <- test_poly_mle %>% dplyr::filter(lambda_r_diff < 10000, lambda_r_diff>-10000, mu_r_diff<10000,mu_r_diff>-10000, cap_r_diff<10000, cap_r_diff>-10000)%>%
  group_by(lambda, .add = T) %>%
  group_by(mu, .add = T) %>%
  group_by(cap, .add = T) %>%
  summarise(mae_lambda = mean(abs(lambda_a_diff)),
            mae_mu = mean(abs(mu_a_diff)),
            mae_cap = mean(abs(cap_a_diff)),
            mean_node = mean(num_nodes)) %>%
  mutate(cap = if_else(cap == 700, true = 800, false = cap))

saveRDS(test_poly_mle2, "C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\robustness\\robustness_MLE_typ.rds")



fill_pal <- "aurora"
color_pal <- "aurora"

p_lambda <- ggplot(test_poly_gnn2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) +
  facet_wrap(.~cap, nrow=1, labeller = as_labeller(~ paste0("Carrying capacity (K): ",.))) +
  geom_text(aes(lambda, mu, label = round(mae_lambda, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
                                            axis.text.x = element_blank()) +
  labs(x=NULL,y="Error λ")

p_mu <- ggplot(test_poly_gnn2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_mu, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
                                            axis.text.x = element_blank(),
                                            strip.text.x = element_blank()) +
  labs(x=NULL,y="Error μ")

p_cap <- ggplot(test_poly_gnn2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_cap, 0))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0), label = scales::label_number(prefix = "λ: ")) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
                                            strip.text.x = element_blank()) +
  labs(x=NULL,y="Error K")

plstm_lambda <- ggplot(test_poly_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) +
  facet_wrap(.~cap, nrow=1, labeller = as_labeller(~ paste0("Carrying capacity (K): ",.))) +
  geom_text(aes(lambda, mu, label = round(mae_lambda, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  labs(x=NULL,y="Error λ")

plstm_mu <- ggplot(test_poly_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_mu, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error μ")

plstm_cap <- ggplot(test_poly_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_cap, 0))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0), label = scales::label_number(prefix = "λ: ")) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error K")


plstmdnn_lambda <- ggplot(test_poly_lstm_dnn2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) +
  facet_wrap(.~cap, nrow=1, labeller = as_labeller(~ paste0("Carrying capacity (K): ",.))) +
  geom_text(aes(lambda, mu, label = round(mae_lambda, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  labs(x=NULL,y="Error λ")

plstmdnn_mu <- ggplot(test_poly_lstm_dnn2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_mu, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error μ")

plstmdnn_cap <- ggplot(test_poly_lstm_dnn2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_cap, 0))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0), label = scales::label_number(prefix = "λ: ")) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error K")

pmleb_lambda <- ggplot(test_poly_mle_opt2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) +
  facet_wrap(.~cap, nrow=1, labeller = as_labeller(~ paste0("Carrying capacity (K): ",.))) +
  geom_text(aes(lambda, mu, label = round(mae_lambda, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  labs(x=NULL,y="Error λ")

pmleb_mu <- ggplot(test_poly_mle_opt2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_mu, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error μ")

pmleb_cap <- ggplot(test_poly_mle_opt2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_cap, 0))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0), label = scales::label_number(prefix = "λ: ")) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error K")

pmlet_lambda <- ggplot(test_poly_mle2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) +
  facet_wrap(.~cap, nrow=1, labeller = as_labeller(~ paste0("Carrying capacity (K): ",.))) +
  geom_text(aes(lambda, mu, label = round(mae_lambda, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  labs(x=NULL,y="Error λ")

pmlet_mu <- ggplot(test_poly_mle2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_mu, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error μ")

pmlet_cap <- ggplot(test_poly_mle2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_cap, 0))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0), label = scales::label_number(prefix = "λ: ")) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error K")

test_poly_gnn2$methods <- "GNN"
test_poly_gnn2$mean_node <- NULL
test_poly_lstm2$methods <- "Boost BT"
test_poly_lstm_dnn2$methods <- "Boost SS+BT"
test_poly_mle2$methods <- "MLE Naive"
test_poly_mle2$mean_node <- NULL
test_poly_mle_opt2$methods <- "MLE Best"
test_poly_mle_opt2$mean_node <- NULL

robustness_combined_df <- rbind(test_poly_gnn2,test_poly_lstm2,test_poly_lstm_dnn2,test_poly_mle2,test_poly_mle_opt2)
robustness_combined_df$mu <- round(robustness_combined_df$mu, 5)
robustness_combined_df$keff <- (robustness_combined_df$lambda - robustness_combined_df$mu) / robustness_combined_df$cap
robustness_combined_df$lambda <- as.factor(robustness_combined_df$lambda)
robustness_combined_df$mu <- as.factor(robustness_combined_df$mu)
robustness_combined_df$cap <- as.factor(robustness_combined_df$cap)

print(arrange(robustness_combined_df, mae_lambda), n = 50)
print(arrange(robustness_combined_df, mae_mu), n = 50)
print(arrange(robustness_combined_df, mae_cap), n = 50)

robustness_combined_df_bak <- robustness_combined_df

##mark_lambda_df <- arrange(robustness_combined_df, mae_lambda)[1,]
##mark_mu_df <- arrange(robustness_combined_df, mae_mu)[24,]
##mark_cap_df <- arrange(robustness_combined_df, mae_cap)[35,]

mark_lambda_df <- robustness_combined_df %>% filter(lambda==1,mu==0.2,cap==800) %>% 
  mutate(text = if_else(methods=="GNN","▼",""),
         type = "NN")
mark_lambda1_df <- rbind(mark_lambda_df, mark_lambda_df %>% mutate(type = "MLE", text = ""))
mark_lambda2_df <- rbind(mark_lambda_df %>% mutate(text = ""),mark_lambda_df %>% mutate(type = "MLE"))

mark_mu_df <- robustness_combined_df %>% filter(lambda==1.5,mu==0.2,cap==800) %>% 
  mutate(text = if_else(methods=="Boost BT","▼",""),
         type = "NN")
mark_mu1_df <- rbind(mark_mu_df, mark_mu_df %>% mutate(type = "MLE", text = ""))
mark_mu_df2 <- robustness_combined_df %>% filter(lambda==3,mu==0.2,cap==800) %>% 
  mutate(text = if_else(methods=="MLE Best","▼",""),
         type = "MLE")
mark_mu2_df <- rbind(mark_mu_df %>% mutate(text = ""),mark_mu_df2 %>% mutate(type = "MLE"))
mark_cap_df <- robustness_combined_df %>% filter(lambda==3,mu==0.2,cap==400) %>% 
  mutate(text = if_else(methods=="GNN","▼",""),
         type = "NN")
mark_cap1_df <- rbind(mark_cap_df, mark_cap_df %>% mutate(type = "MLE", text = ""))
mark_cap_df2 <- robustness_combined_df %>% filter(lambda==3,mu==0.2,cap==200) %>% 
  mutate(text = if_else(methods=="MLE Best","▼",""),
         type = "MLE")
mark_cap2_df <- rbind(mark_cap_df %>% mutate(text = ""),mark_cap_df2 %>% mutate(type = "MLE"))


`%||%` <- function(a, b) {
  if(is.null(a)) b else a
}

GeomRtile <- ggproto("GeomRtile", 
                     statebins:::GeomRrect, # 1) only change compared to ggplot2:::GeomTile
                     
                     extra_params = c("na.rm"),
                     setup_data = function(data, params) {
                       data$width <- data$width %||% params$width %||% resolution(data$x, FALSE)
                       data$height <- data$height %||% params$height %||% resolution(data$y, FALSE)
                       
                       transform(data,
                                 xmin = x - width / 2,  xmax = x + width / 2,  width = NULL,
                                 ymin = y - height / 2, ymax = y + height / 2, height = NULL
                       )
                     },
                     default_aes = aes(
                       fill = "grey20", colour = NA, size = 0.1, linetype = 1,
                       alpha = NA, width = NA, height = NA
                     ),
                     required_aes = c("x", "y"),
                     
                     # These aes columns are created by setup_data(). They need to be listed here so
                     # that GeomRect$handle_na() properly removes any bars that fall outside the defined
                     # limits, not just those for which x and y are outside the limits
                     non_missing_aes = c("xmin", "xmax", "ymin", "ymax"),
                     draw_key = draw_key_polygon
)

geom_rtile <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       radius = grid::unit(6, "pt"), # 2) add radius argument
                       ...,
                       linejoin = "mitre",
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRtile, # 3) use ggproto object here
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = rlang::list2(
      radius = radius,
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

oldK <- GeomText$draw_key # to save for later

# define new key
# if you manually add colours then add vector of colours 
# instead of `scales::hue_pal()(length(var))`
GeomText$draw_key <- function (data, params, size, 
                               var=c("▼","▼"), 
                               cols=c("blue","salmon")) {
  
  # sort as ggplot sorts these alphanumerically / or levels of factor
  txt <- if(is.factor(var)) levels(var) else sort(var)
  txt <- txt[match(data$colour, cols)]
  
  textGrob(txt, 0.5, 0.5,  
           just="center", 
           gp = gpar(col = alpha(data$colour, data$alpha), 
                     fontfamily = data$family, 
                     fontface = data$fontface, 
                     fontsize = data$size * .pt))
}

# reset key
GeomText$draw_key <- oldK


fill_pal <- "lake_superior"
fill_pal <- "lumina"
fill_pal <- "frost"
fill_pal <- "afternoon_prarie"
fill_pal2 <- "aurora" 
fill_pal2 <- "afternoon_prarie" 

cust_labeller_lambda <- function(x) paste0("λ: ", x)
cust_labeller_cap <- function(x) paste0("K: ", x)

ddd_new_robustness_lambda <- ggplot(robustness_combined_df) + geom_rtile(aes(mu, mae_lambda, fill=keff),color = "#666666",lty="dashed",
                                                                      height = Inf,col=NA) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K", reverse = T) +
  new_scale_fill() +
  geom_bar(aes(mu, mae_lambda, fill = methods, group= methods), color="#666666",
           stat = "identity", position = "dodge") + 
  facet_grid(lambda~cap, labeller = labeller(
    lambda = cust_labeller_lambda,
    cap = cust_labeller_cap
  )) +
  scale_x_discrete(expand=c(0,0),label = function(x) paste0("μ: ",x))+
  scale_y_continuous(trans = "sqrt", breaks = c(0,0.1,0.5,2.0))+
  scale_fill_nord(fill_pal2,discrete = T, name="Method") +
  theme(axis.line = element_blank(), axis.ticks.x = element_blank(),
        plot.background = element_blank(), panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x=NULL,y="Mean absolute error (speciation rate λ)") +
  geom_hline(yintercept=0.5,color="black",lty="twodash",size=0.1)+
  new_scale_color()+
  scale_color_manual(values = c("NN"="salmon","MLE"="blue"), labels = c("MLE"="Best of all", "NN"="Best of NN")) + 
  geom_text(data = mark_lambda1_df, aes(mu, mae_lambda, group= methods, label=text, color = type), 
            size = 5, vjust = -0.15, position = position_dodge(width=.9))+
  geom_text(data = mark_lambda2_df, aes(mu, mae_lambda, group= methods, label=text, color = type), 
            size = 5, vjust = -1.2, position = position_dodge(width=.9), show.legend = F)+
  labs(color="Mark")


ddd_new_robustness_mu <- ggplot(robustness_combined_df) + geom_rtile(aes(mu, mae_mu, fill=keff),color = "#666666",lty="dashed",
                                                                      height = Inf,col=NA) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K", reverse = T) +
  new_scale_fill() +
  geom_bar(aes(mu, mae_mu, fill = methods, group= methods), color="#666666",
           stat = "identity", position = "dodge") + 
  facet_grid(lambda~cap, labeller = labeller(
    lambda = cust_labeller_lambda,
    cap = cust_labeller_cap
  )) +
  scale_x_discrete(expand=c(0,0),label = function(x) paste0("μ: ",x))+
  scale_y_continuous(trans = "sqrt", breaks = c(0,0.03,0.3,1.0))+
  scale_fill_nord(fill_pal2,discrete = T, name="Method") +
  theme(axis.line = element_blank(), axis.ticks.x = element_blank(),
        plot.background = element_blank(), panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x=NULL,y="Mean absolute error (extinction rate μ)") +
  geom_hline(yintercept=0.3,color="black",lty="dotted",size=0.1)+
  new_scale_color()+
  scale_color_manual(values = c("NN"="salmon","MLE"="blue"), labels = c("MLE"="Best of all", "NN"="Best of NN")) + 
  geom_text(data = mark_mu1_df, aes(mu, mae_mu, group= methods, label=text, color = type), 
            size = 5, vjust = -0.15, position = position_dodge(width=.9))+
  geom_text(data = mark_mu2_df, aes(mu, mae_mu, group= methods, label=text, color = type), 
            size = 5, vjust = -0.15, position = position_dodge(width=.9), show.legend = F)+
  labs(color="Mark")

ddd_new_robustness_cap <- ggplot(robustness_combined_df) + geom_rtile(aes(mu, mae_cap, fill=keff),color = "#666666",lty="dashed",
                                           height = Inf,col=NA) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K", reverse = T) +
  new_scale_fill() +
  geom_bar(aes(mu, mae_cap, fill = methods, group= methods), color="#666666",
                                          stat = "identity", position = "dodge") + 
  facet_grid(lambda~cap, labeller = labeller(
    lambda = cust_labeller_lambda,
    cap = cust_labeller_cap
  )) +
  scale_x_discrete(expand=c(0,0),label = function(x) paste0("μ: ",x))+
  scale_y_continuous(trans = "sqrt", breaks = c(0,20,100,500, 2000))+
  scale_fill_nord(fill_pal2,discrete = T, name="Method") +
  theme(axis.line = element_blank(), axis.ticks.x = element_blank(),
        plot.background = element_blank(), panel.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x=NULL,y="Mean absolute error (carrying capacity K)") +
  geom_hline(yintercept=500,color="black",lty="twodash",size=0.1)+
  new_scale_color()+
  scale_color_manual(values = c("NN"="salmon","MLE"="blue"), labels = c("MLE"="Best of all", "NN"="Best of NN")) + 
  geom_text(data = mark_cap1_df, aes(mu, mae_cap, group= methods, label=text, color = type), 
            size = 5, vjust = -0.15, position = position_dodge(width=.9))+
  geom_text(data = mark_cap2_df, aes(mu, mae_cap, group= methods, label=text, color = type), 
            size = 5, vjust = -0.15, position = position_dodge(width=.9), show.legend = F)+
  labs(color="Mark")

ddd_new_robustness <- ddd_new_robustness_lambda + ddd_new_robustness_mu + ddd_new_robustness_cap + 
  plot_layout(ncol = 1, guides = "collect") + plot_annotation(title = "Robustness Analysis DDD between Estimation Methods against Carrying Capacity Effect Strength")

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_new_robustness2.png",
       ddd_new_robustness,
       width = 16,
       height = 19,
       dpi = 300,
       device = "png")

ddd_robustness <- grid.arrange(
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Neural Network (GNN) Robustness") +
      inset_element(
        p_lambda + p_mu + p_cap +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Neural Network (Boosting BT) Robustness") +
      inset_element(
        plstm_lambda + plstm_mu + plstm_cap +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Neural Network (Boosting SS + BT) Robustness") +
      inset_element(
        plstmdnn_lambda + plstmdnn_mu + plstmdnn_cap +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Maximum Likelihood Estimation (Naive Case) Robustness") +
      inset_element(
        pmlet_lambda + pmlet_mu + pmlet_cap +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Maximum Likelihood Estimation (Best Case) Robustness") +
      inset_element(
        pmleb_lambda + pmleb_mu + pmleb_cap +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  ncol = 1
)

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_robustness.png",
       ddd_robustness,
       width = 12,
       height = 16,
       dpi = "retina",
       device = "png")

p_lambda_dnn <- ggplot(test_poly_dnn_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) +
  facet_wrap(.~cap, nrow=1, labeller = as_labeller(~ paste0("Carrying capacity (K): ",.))) +
  geom_text(aes(lambda, mu, label = round(mae_lambda_dnn, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  labs(x=NULL,y="Error λ")

p_mu_dnn <- ggplot(test_poly_dnn_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_mu_dnn, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error μ")

p_cap_dnn <- ggplot(test_poly_dnn_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_cap_dnn, 0))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0), label = scales::label_number(prefix = "λ: ")) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error K")

p_lambda_lstm <- ggplot(test_poly_dnn_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) +
  facet_wrap(.~cap, nrow=1, labeller = as_labeller(~ paste0("Carrying capacity (K): ",.))) +
  geom_text(aes(lambda, mu, label = round(mae_lambda_lstm, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  labs(x=NULL,y="Error λ")

p_mu_lstm <- ggplot(test_poly_dnn_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_mu_lstm, 2))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error μ")

p_cap_lstm <- ggplot(test_poly_dnn_lstm2) + geom_tile(aes(lambda, mu, fill = (lambda-mu)/cap)) + facet_wrap(.~cap, nrow=1) +
  geom_text(aes(lambda, mu, label = round(mae_cap_lstm, 0))) +
  scale_fill_nord(fill_pal,discrete = F, trans="log10", name="(λ - μ) / K") +
  scale_color_nord(color_pal,discrete = F) +
  scale_x_continuous(expand=c(0,0), label = scales::label_number(prefix = "λ: ")) +
  scale_y_continuous(expand=c(0,0), breaks = c(0.2, 0.4, 0.6, 0.8), label = scales::label_number(prefix = "μ: ")) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        strip.text.x = element_blank()) +
  labs(x=NULL,y="Error K")

ddd_robustness_dnn_lstm <- grid.arrange(
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Neural Network (GNN) Robustness") +
      inset_element(
        p_lambda + p_mu + p_cap +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Neural Network (DNN) Robustness") +
      inset_element(
        p_lambda_dnn + p_mu_dnn + p_cap_dnn +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  patchworkGrob(
    ggplot() +
      annotation_custom(grid::roundrectGrob(
        r = unit(0.03, "snpc"),
        gp = grid::gpar(lty = "dashed", lwd = 1, fill = NA, col = "gray")
      )) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(plot.title = element_text(margin = margin(b = 5.5))) +
      labs(title = "Neural Network (LSTM) Robustness") +
      inset_element(
        p_lambda_lstm + p_mu_lstm + p_cap_lstm +
          plot_layout(ncol = 1, guides = "collect") +
          plot_annotation(title = "Neural Network"),
        left = 0,
        right = 1,
        top = 1,
        bottom = 0,
        on_top = FALSE
      )
  ),
  ncol = 1
)

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\ddd_robustness_dnn_lstm.png",
       ddd_robustness_dnn_lstm,
       width = 12,
       height = 10,
       dpi = "retina",
       device = "png")

test <- readRDS("D:\\Data\\Empirical Trees\\Export\\GNN\\tree\\ST\\st_Amphibia_Alsodidae.rds")


condamine_mle <- eveGNN::load_empirical_mle_result("D:\\Habrok\\Data\\EMP")

condamine_mle_ddd <- condamine_mle$DDD

condamine_gnn_ddd <- readRDS("D:\\Habrok\\Data\\EMP\\empirical_gnn_2_lstm_result.rds")

clean_text <- function(text) {
  gsub("\\['(.*?)'\\]", "\\1", text)
}

condamine_gnn_ddd$family <- sapply(condamine_gnn_ddd$family, clean_text)
condamine_gnn_ddd$tree <- sapply(condamine_gnn_ddd$tree, clean_text)

condamine_gnn_ddd <- condamine_gnn_ddd %>% rename(Tree = tree, Family = family)

condamine_data <- left_join(condamine_gnn_ddd, condamine_mle_ddd, by=c("Family", "Tree"))

condamine_data$num_nodes <- as.numeric(condamine_data$num_nodes)
condamine_data$pred_cap_before <- condamine_data$pred_cap_before * 1000
condamine_data$pred_cap_after <- condamine_data$pred_cap_after * 1000

condamine_data2 <- condamine_data %>% mutate(net_div_nn=pred_lambda_after - pred_mu_after,
                                             net_div_mle=lambda - mu,
                                             sp_nn=pred_lambda_after,
                                             ext_nn=pred_mu_after,
                                             cc_nn=pred_cap_after,
                                             sp_mle=lambda,
                                             ext_mle=mu,
                                             cc_mle=cap,
                                             cef_nn=(pred_lambda_after-pred_mu_after)/pred_cap_after,
                                             cef_mle=(lambda-mu)/cap) %>%
  select(-pred_lambda_after,-pred_lambda_before,-pred_mu_after,-pred_mu_before,-pred_cap_after,-pred_cap_before,-lambda,-mu,-cap,-df,-conv,-Model)


ggthemr::ggthemr("flat")

color_palette_emp <- "lumina"
range_color_scale <- range(condamine_data2$num_nodes)

pemp_nd_ddd <- ggplot(condamine_data2) + ggpointdensity::geom_pointdensity(aes(x = net_div_nn, y = net_div_mle, color = num_nodes)) +
  lims(x=c(0,3),y=c(0,3)) + geom_abline(slope=1,intercept=0,linewidth=1,linetype="dashed",color="grey") +
  labs(x="Boosting BT", y = "MLE", title = "Net diversification rate (λ - μ)") +
  scale_color_nord(palette = color_palette_emp,discrete=FALSE, limits=range_color_scale, breaks=c(200,500,1000), name="Size") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1/1)

pemp_sp_ddd <- ggplot(condamine_data2) + ggpointdensity::geom_pointdensity(aes(x = sp_nn, y = sp_mle, color = num_nodes)) +
  lims(x=c(0,3),y=c(0,3)) + geom_abline(slope=1,intercept=0,linewidth=1,linetype="dashed",color="grey") +
  labs(x="Boosting BT", y = "MLE", title = "Speciation rate (λ)") +
  scale_color_nord(palette = color_palette_emp,discrete=FALSE, limits=range_color_scale, breaks=c(200,500,1000), name="Size") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1/1)

pemp_ext_ddd <- ggplot(condamine_data2) + ggpointdensity::geom_pointdensity(aes(x = ext_nn, y = ext_mle, color = num_nodes)) +
  lims(x=c(0,3),y=c(0,3)) + geom_abline(slope=1,intercept=0,linewidth=1,linetype="dashed",color="grey") +
  labs(x="Boosting BT", y = "MLE", title = "Extinction rate (μ)") +
  scale_color_nord(palette = color_palette_emp,discrete=FALSE, limits=range_color_scale, breaks=c(200,500,1000), name="Size") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1/1)

pemp_cc_ddd <- ggplot(condamine_data2) + ggpointdensity::geom_pointdensity(aes(x = cc_nn, y = cc_mle, color = num_nodes)) +
  lims(x=c(0,1000),y=c(0,1000)) + geom_abline(slope=1,intercept=0,linewidth=1,linetype="dashed",color="grey") +
  labs(x="Boosting BT", y = "MLE", title = "Carrying capacity (K)") +
  scale_color_nord(palette = color_palette_emp,discrete=FALSE, limits=range_color_scale, breaks=c(200,500,1000), name="Size") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1/1)

emp_est_plot <- pemp_nd_ddd+pemp_sp_ddd+pemp_ext_ddd+pemp_cc_ddd+plot_layout(nrow = 1,guides = "collect") +
  plot_annotation(title="Parameter Estimation on Empirical Trees (Diversity-Dependent Diversification Scenario), Neural Network against MLE")

ggsave("ddd_empirical.png",
       emp_est_plot,
       width=20,
       height=5,
       dpi="retina",
       device="png")

# ▒▒▒ 0. set theme and palette --------------------------------------------
ggthemr::ggthemr("flat")
color_palette_emp <- "lumina"
range_color_scale <- range(condamine_data2$num_nodes)

# ▒▒▒ 1. helper for scatter‑plots coloured by size -------------------------
scatter_fun <- function(df, x, y, ttl, x_lims, y_lims) {
  ggplot(df) +
    ggpointdensity::geom_pointdensity(
      aes(x = .data[[x]], y = .data[[y]], colour = num_nodes)
    ) +
    lims(x = x_lims, y = y_lims) +
    geom_abline(slope = 1, intercept = 0, linewidth = 1,
                linetype = "dashed", colour = "grey") +
    labs(x = "Boosting BT", y = "MLE", title = ttl) +
    scale_color_nord(palette = color_palette_emp, discrete = FALSE,
                     limits = range_color_scale,
                     breaks  = c(200, 500, 1000),
                     name    = "Size") +
    theme(plot.background  = element_blank(),
          panel.background = element_blank(),
          aspect.ratio     = 1)
}

# ▒▒▒ 2. helper for |Δ| vs size panels ------------------------------------
err_fun <- function(df, nn_col, mle_col, ttl, ylim) {
  ggplot(df, aes(x = num_nodes,
                 y = abs(.data[[nn_col]] - .data[[mle_col]]))) +
    geom_point(alpha = .35, size = 1) +
    #geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
    scale_y_continuous(limits = ylim, name = "|Δ|") +
    scale_x_continuous(name = "Tree size") +
    labs(title = ttl) +
    theme_minimal(base_size = 9) +
    theme(plot.title      = element_text(hjust = .5),
          legend.position = "none",
          aspect.ratio    = .5)
}

# ▒▒▒ 3. build the eight panels -------------------------------------------
p_nd  <- scatter_fun(condamine_data2, "net_div_nn", "net_div_mle",
                     "Net diversification rate (λ − μ)",
                     c(0, 3), c(0, 3))
p_sp  <- scatter_fun(condamine_data2, "sp_nn", "sp_mle",
                     "Speciation rate (λ)", c(0, 3), c(0, 3))
p_ext <- scatter_fun(condamine_data2, "ext_nn", "ext_mle",
                     "Extinction rate (μ)", c(0, 3), c(0, 3))
p_cc  <- scatter_fun(condamine_data2, "cc_nn", "cc_mle",
                     "Carrying capacity (K)", c(0, 1000), c(0, 1000))

e_nd  <- err_fun(condamine_data2, "net_div_nn", "net_div_mle",
                 "|λ − μ| vs size", c(0, 3))
e_sp  <- err_fun(condamine_data2, "sp_nn", "sp_mle",
                 "|λ| vs size", c(0, 3))
e_ext <- err_fun(condamine_data2, "ext_nn", "ext_mle",
                 "|μ| vs size", c(0, 3))
e_cc  <- err_fun(condamine_data2, "cc_nn", "cc_mle",
                 "|K| vs size", c(0, 1000))

# ▒▒▒ 4. assemble ----------------------------------------------------------
emp_est_plot <- (p_nd | p_sp | p_ext | p_cc) /
  (e_nd | e_sp | e_ext | e_cc) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Parameter Estimation on Empirical Trees (Diversity-Dependent Diversification Scenario), Neural Network against MLE"
  )

# ▒▒▒ 5. export ------------------------------------------------------------
ggsave("ddd_empirical_size_effect.png",
       emp_est_plot,
       width  = 20,
       height = 10,   # doubled to accommodate 2 rows
       dpi    = "retina")


## PBD Boosting Results
pbd_boost_dnn <- readRDS("D:/Habrok/Data/STBO/PBD_FREE_TES/PBD_FREE_TES_gnn_2_dnn_compensation.rds")
pbd_boost_dnn <- pbd_boost_dnn %>% dplyr::mutate_all(as.numeric)
pbd_boost_dnn <- pbd_boost_dnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                          num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                          num_nodes >= 500 ~ "Large"))
pbd_boost_dnn_sampled <- pbd_boost_dnn[sample(1:nrow(pbd_boost_dnn), 2000),]

pbd_boost_dnn_sampled <- pbd_boost_dnn_sampled %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(pred_lambda3_after,
                                                       pred_lambda2_after,
                                                       pred_mu2_after))) %>%
  dplyr::ungroup()

pbd_boost_lstm <- readRDS("D:/Habrok/Data/STBO/PBD_FREE_TES/PBD_FREE_TES_gnn_2_lstm_compensation.rds")
pbd_boost_lstm <- pbd_boost_lstm %>% dplyr::mutate_all(as.numeric)
pbd_boost_lstm <- pbd_boost_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                            num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                            num_nodes >= 500 ~ "Large"))
pbd_boost_lstm_sampled <- pbd_boost_lstm[sample(1:nrow(pbd_boost_lstm), 2000),]

pbd_boost_lstm_sampled <- pbd_boost_lstm_sampled %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(pred_lambda3_after,
                                                       pred_lambda2_after,
                                                       pred_mu2_after))) %>%
  dplyr::ungroup()

pbd_boost_lstm_after_dnn <- readRDS("D:/Habrok/Data/STBO/PBD_FREE_TES/PBD_FREE_TES_gnn_2_lstm_compensation_after_dnn.rds")
pbd_boost_lstm_after_dnn <- pbd_boost_lstm_after_dnn %>% dplyr::mutate_all(as.numeric)
pbd_boost_lstm_after_dnn <- pbd_boost_lstm_after_dnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                                num_nodes >= 500 ~ "Large"))
pbd_boost_lstm_after_dnn_sampled <- pbd_boost_lstm_after_dnn[sample(1:nrow(pbd_boost_lstm_after_dnn), 2000),]

pbd_boost_lstm_after_dnn_sampled <- pbd_boost_lstm_after_dnn_sampled %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(pred_lambda3_after_lstm_after_dnn,
                                                       pred_lambda2_after_lstm_after_dnn,
                                                       pred_mu2_after_lstm_after_dnn))) %>%
  dplyr::ungroup()

pbd_boost_dnn_after_lstm <- readRDS("D:/Habrok/Data/STBO/PBD_FREE_TES/PBD_FREE_TES_gnn_2_dnn_compensation_after_lstm.rds")
pbd_boost_dnn_after_lstm <- pbd_boost_dnn_after_lstm %>% dplyr::mutate_all(as.numeric)
pbd_boost_dnn_after_lstm <- pbd_boost_dnn_after_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                                num_nodes >= 500 ~ "Large"))
pbd_boost_dnn_after_lstm_sampled <- pbd_boost_dnn_after_lstm[sample(1:nrow(pbd_boost_dnn_after_lstm), 2000),]

pbd_boost_dnn_after_lstm_sampled <- pbd_boost_dnn_after_lstm_sampled %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(pred_lambda3_after_dnn_after_lstm,
                                                       pred_lambda2_after_dnn_after_lstm,
                                                       pred_mu2_after_dnn_after_lstm))) %>%
  dplyr::ungroup()

pbd_boost_dnn_and_lstm <- readRDS("D:/Habrok/Data/STBO/PBD_FREE_TES/PBD_FREE_TES_gnn_2_lstm_dnn_compensation.rds")
pbd_boost_dnn_and_lstm <- pbd_boost_dnn_and_lstm %>% dplyr::mutate_all(as.numeric)
pbd_boost_dnn_and_lstm <- pbd_boost_dnn_and_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                                num_nodes >= 500 ~ "Large"))
pbd_boost_dnn_and_lstm_sampled <- pbd_boost_dnn_and_lstm[sample(seq_len(nrow(pbd_boost_dnn_and_lstm)), 2000),]

pbd_boost_dnn_and_lstm_sampled <- pbd_boost_dnn_and_lstm_sampled %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(pred_lambda3_after_lstm_after_dnn,
                                                       pred_lambda2_after_lstm_after_dnn,
                                                       pred_mu2_after_lstm_after_dnn))) %>%
  dplyr::ungroup()

pbd_boost_gnn <- load_final_difference_by_layer(path = "D:/Habrok/Data/STBO",
                                                task_type = "PBD_FREE_TES", model = "diffpool", depth = 2)
pbd_boost_gnn <- pbd_boost_gnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                          num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                          num_nodes >= 500 ~ "Large"))
pbd_boost_gnn_sampled <- pbd_boost_gnn[sample(1:nrow(pbd_boost_gnn), 2000),]

pbd_boost_gnn_sampled <- pbd_boost_gnn_sampled %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(lambda3_pred, lambda2_pred, mu2_pred))) %>%
  dplyr::ungroup()

pbd_mle_opt <- load_separated_mle_result(path = "D:/Habrok/Data/MLE", task_type = "PBD", model_type = "diffpool", no_init = FALSE)

pbd_mle_opt <- pbd_mle_opt %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                      num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                      num_nodes >= 500 ~ "Large"))

pbd_mle_opt <- pbd_mle_opt %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(lambda3_pred, lambda2_pred, mu2_pred))) %>%
  dplyr::ungroup()

pbd_mle_typ <- load_separated_mle_result(path = "D:/Habrok/Data/MLE", task_type = "PBD", model_type = "diffpool", no_init = TRUE)

pbd_mle_typ <- pbd_mle_typ %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                      num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                      num_nodes >= 500 ~ "Large"))

pbd_mle_typ <- pbd_mle_typ %>%
  dplyr::rowwise() %>%
  dplyr::mutate(durspec = PBD::pbd_durspec_mean(c(lambda3, lambda2, mu2)),
                durspec_pred = PBD::pbd_durspec_mean(c(lambda3_pred, lambda2_pred, mu2_pred))) %>%
  dplyr::ungroup()

point_alpha <- 0.4
df.borders <- data.frame(intercept=c(0,0,10000),slope=c(1,0,0),Reference=c('zero','same', "mean"))

plot_list_pbd_true <- list()
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda1, lambda1 - lambda1_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="GNN      \nError", title = bquote(italic(lambda)[1]))
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda2, lambda2 - lambda2_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2)+ ggplot2::ylim(-15,15) + ggplot2::labs(x = "True parameter value", y = NULL, title = bquote(italic(lambda)[2]))
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda3, lambda3 - lambda3_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL, title = bquote(italic(lambda)[3]))
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu1, mu1 - mu1_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL, title = bquote(italic(mu)[1])) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu2, mu2 - mu2_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL, title = bquote(italic(mu)[2])) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,3) + ggplot2::labs(x = "True parameter value", y = NULL, title = bquote(italic(tau))) + ggplot2::xlim(0,3)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,40) + ggplot2::labs(x = "True parameter value", y = NULL, title = bquote(italic(tau))) + ggplot2::xlim(3,40)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda1, lambda1 - pred_lambda1_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="Boost SS      \nError")
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda2, lambda2 - pred_lambda2_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2)+ ggplot2::ylim(-15,15) + ggplot2::labs(x = "True parameter value", y = NULL)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda3, lambda3 - pred_lambda3_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu1, mu1 - pred_mu1_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu2, mu2 - pred_mu2_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,3) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,3)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,40) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(3,40)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda1, lambda1 - pred_lambda1_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="Boost BT      \nError")
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda2, lambda2 - pred_lambda2_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2)+ ggplot2::ylim(-15,15) + ggplot2::labs(x = "True parameter value", y = NULL)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda3, lambda3 - pred_lambda3_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(mu1, mu1 - pred_mu1_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(mu2, mu2 - pred_mu2_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,3) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,3)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,40) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(3,40)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda1, lambda1 - pred_lambda1_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="Boost SS+BT      \nError")
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda2, lambda2 - pred_lambda2_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2)+ ggplot2::ylim(-15,15) + ggplot2::labs(x = "True parameter value", y = NULL)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda3, lambda3 - pred_lambda3_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu1, mu1 - pred_mu1_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu2, mu2 - pred_mu2_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.64 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,3) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,3)
plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2)+ ggplot2::ylim(-1,40) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(3,40)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(lambda1, lambda1_a_diff, color = small), alpha = point_alpha) +ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="MLE Naive      \nError")

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(lambda2, lambda2_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-15,15) + ggplot2::labs(x = "True parameter value", y = NULL)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(lambda3, lambda3_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(mu1, mu1_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(mu2, mu2_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,3) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,3)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,40) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(3,40)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(lambda1, lambda1_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="MLE Best      \nError")

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(lambda2, lambda2_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-15,15) + ggplot2::labs(x = "True parameter value", y = NULL)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(lambda3, lambda3_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::theme(axis.title.x = element_text(hjust = 0.5))

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(mu1, mu1_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(mu2, mu2_a_diff, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.64)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,3) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,3)

plot_list_pbd_true[[length(plot_list_pbd_true) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, durspec - durspec_pred, color = small), alpha = point_alpha) + ggplot2::labs(color="Tree size") +ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+ggplot2::ylim(-1,40) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(3,40)


for (i in 1:(length(plot_list_pbd_true)-14)) {
  plot_list_pbd_true[[i]] <- plot_list_pbd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list_pbd_true)-13):length(plot_list_pbd_true)) {
  plot_list_pbd_true[[i]] <- plot_list_pbd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

pbd_plot_true_value <- patchwork::wrap_plots(plot_list_pbd_true) +
  patchwork::plot_layout(ncol = 7, guides = "collect", axis_title = "collect_x") + patchwork::plot_annotation(title = "Performance Analysis PBD against True Value")

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\pbd_plot_true_value2.png",
       pbd_plot_true_value,
       width = 20,
       height = 10,
       dpi = "retina",
       device = "png")
ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\pbd_plot_true_value2.png",
       pbd_plot_true_value,
       width = 20,
       height = 10,
       dpi = "retina",
       device = "png")

## PBD durspec vs treesize
pbd_durspec_size <- ggplot2::ggplot(pbd_boost_gnn_sampled) +
  ggplot2::geom_point(ggplot2::aes(durspec, num_nodes, color = small)) +
  ggplot2::theme(panel.background = ggplot2::element_blank()) +
  ggplot2::labs(x="Mean duration of speciation", y="Tree size", color="Size")

pbd_durspec_lambda1 <- ggplot2::ggplot(pbd_boost_gnn_sampled) +
  ggplot2::geom_point(ggplot2::aes(durspec, lambda1, color = small)) +
  ggplot2::theme(panel.background = ggplot2::element_blank()) +
  ggplot2::labs(x="Mean duration of speciation", y="Speciation good species", color="Size")

pbd_durspec_lambda2 <- ggplot2::ggplot(pbd_boost_gnn_sampled) +
  ggplot2::geom_point(ggplot2::aes(durspec, lambda2, color = small)) +
  ggplot2::theme(panel.background = ggplot2::element_blank()) +
  ggplot2::labs(x="Mean duration of speciation", y="Speciation completion rate", color="Size")

pbd_durspec_lambda3 <- ggplot2::ggplot(pbd_boost_gnn_sampled) +
  ggplot2::geom_point(ggplot2::aes(durspec, lambda3, color = small)) +
  ggplot2::theme(panel.background = ggplot2::element_blank()) +
  ggplot2::labs(x="Mean duration of speciation", y="Speciation incipient species", color="Size")

pbd_durspec_mu1 <- ggplot2::ggplot(pbd_boost_gnn_sampled) +
  ggplot2::geom_point(ggplot2::aes(durspec, mu1, color = small)) +
  ggplot2::theme(panel.background = ggplot2::element_blank()) +
  ggplot2::labs(x="Mean duration of speciation", y="Extinction good species", color="Size")

pbd_durspec_mu2 <- ggplot2::ggplot(pbd_boost_gnn_sampled) +
  ggplot2::geom_point(ggplot2::aes(durspec, mu2, color = small)) +
  ggplot2::theme(panel.background = ggplot2::element_blank()) +
  ggplot2::labs(x="Mean duration of speciation", y="Extinction incipient species", color="Size")

pbd_durspec_pars <- pbd_durspec_size + pbd_durspec_lambda1 + pbd_durspec_lambda2 + pbd_durspec_lambda3 + pbd_durspec_mu1 + pbd_durspec_mu2 +
  patchwork::plot_layout(ncol = 3, guides = "collect", axis_titles = "collect") +
  patchwork::plot_annotation(title = "Correlation between the Mean Duration of Speciation and the Parameters")

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\pbd_durspec_pars.png",
       pbd_durspec_pars,
       width = 12,
       height = 6,
       dpi = "retina",
       device = "png")


## PBD results tree size
nord_palette <-  "lumina"
mae_digits_float <- 2
mae_digits_int <- 0
mae_text_size <- 5
mae_vjust <- 1.6
mae_color <- "#666666"
x_breaks <- c(0, 200, 500, 1000, 1500, 2000)
range_color_scale <- c(0,40)

ggthemr::ggthemr("flat")
plot_list_pbd_size <- list()
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda1 - lambda1_pred, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y="GNN      \nError", title = bquote(italic(lambda)[1]))
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda2 - lambda2_pred, color = durspec)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(lambda)[2]))
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda3 - lambda3_pred, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(lambda)[3]))
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu1 - mu1_pred, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(mu)[1]))
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu2 - mu2_pred, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(mu)[2]))
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda1 - pred_lambda1_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y="Boost SS      \nError")
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda2 - pred_lambda2_after, color = durspec)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda3 - pred_lambda3_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu1 - pred_mu1_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu2 - pred_mu2_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda1 - pred_lambda1_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y="Boost BT      \nError")
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda2 - pred_lambda2_after, color = durspec)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda3 - pred_lambda3_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu1 - pred_mu1_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu2 - pred_mu2_after, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda1 - pred_lambda1_after_lstm_after_dnn, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y="Boost SS+BT      \nError")
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda2 - pred_lambda2_after_lstm_after_dnn, color = durspec)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda3 - pred_lambda3_after_lstm_after_dnn, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu1 - pred_mu1_after_lstm_after_dnn, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu2 - pred_mu2_after_lstm_after_dnn, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda1_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y="MLE Naive      \nError")
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda2_a_diff, color = durspec)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda3_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu1_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu2_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda1_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y="MLE Best      \nError")
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda2_a_diff, color = durspec)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda3_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes",y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu1_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_size[[length(plot_list_pbd_size) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu2_a_diff, color = durspec)) + ggplot2::ylim(-2,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL)

for (i in 1:length(plot_list_pbd_size)) {
  plot_list_pbd_size[[i]] <- plot_list_pbd_size[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 200, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 500, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 2000, color = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = x_breaks, lim = c(0,2100)) +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))+
    ggplot2::labs(color = "τ")
}

pbd_plot_tree_size <- patchwork::wrap_plots(plot_list_pbd_size) +
  patchwork::plot_layout(ncol = 5, guides = "collect") +
  patchwork::plot_annotation(title = "Performance Analysis PBD against Phylogeny Size (Absolute Error)")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\pbd_plot_tree_size2.png",
                pbd_plot_tree_size,
                device = "png", width = 16,height = 9,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\pbd_plot_tree_size2.png",
                pbd_plot_tree_size,
                device = "png", width = 16,height = 9,dpi = "retina")

## PBD durspec
range_color_scale2 <- c(0,2100)
x_breaks2 <- c(0, 10, 20, 30, 40)
ggthemr::ggthemr("flat")
plot_list_pbd_durspec <- list()
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1 - lambda1_pred, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="GNN      \nError", title = bquote(italic(lambda)[1]))
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2 - lambda2_pred, color = num_nodes)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(lambda)[2]))
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3 - lambda3_pred, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(lambda)[3]))
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu1 - mu1_pred, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(mu)[1]))
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu2 - mu2_pred, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(mu)[2]))
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1 - pred_lambda1_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="Boost SS      \nError")
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2 - pred_lambda2_after, color = num_nodes)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3 - pred_lambda3_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu1 - pred_mu1_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu2 - pred_mu2_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1 - pred_lambda1_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="Boost BT      \nError")
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2 - pred_lambda2_after, color = num_nodes)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3 - pred_lambda3_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu1 - pred_mu1_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu2 - pred_mu2_after, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1 - pred_lambda1_after_lstm_after_dnn, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="Boost SS+BT      \nError")
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2 - pred_lambda2_after_lstm_after_dnn, color = num_nodes)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3 - pred_lambda3_after_lstm_after_dnn, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu1 - pred_mu1_after_lstm_after_dnn, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, mu2 - pred_mu2_after_lstm_after_dnn, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="MLE Naive      \nError")
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2_a_diff, color = num_nodes)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, mu1_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, mu2_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="MLE Best      \nError")
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2_a_diff, color = num_nodes)) + ggplot2::ylim(-30,15) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = "Mean duration of speciation",y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, mu1_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec[[length(plot_list_pbd_durspec) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, mu2_a_diff, color = num_nodes)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)

for (i in 1:length(plot_list_pbd_durspec)) {
  plot_list_pbd_durspec[[i]] <- plot_list_pbd_durspec[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))+
    ggplot2::scale_x_continuous(breaks = x_breaks2, lim = c(0,40)) +
    ggplot2::labs(color = "Tree size")
}

pbd_plot_durspec <- patchwork::wrap_plots(plot_list_pbd_durspec) +
  patchwork::plot_layout(ncol = 5, guides = "collect") +
  patchwork::plot_annotation(title = "Performance Analysis PBD against Mean Duration of Speciation (Absolute Error)")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\pbd_plot_durspec.png",
                pbd_plot_durspec,
                device = "png", width = 16,height = 9,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\pbd_plot_durspec2.png",
                pbd_plot_durspec,
                device = "png", width = 16,height = 9,dpi = "retina")

## PBD durspec relative
range_color_scale2 <- c(0,2100)
x_breaks2 <- c(0, 10, 20, 30, 40)
ggthemr::ggthemr("flat")
plot_list_pbd_durspec_rel <- list()
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda1 - lambda1_pred) / lambda1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="GNN      \nError", title = bquote(italic(lambda)[1]))
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda2 - lambda2_pred) / lambda2, color = num_nodes)) + ggplot2::ylim(-1000,5) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(lambda)[2]))
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda3 - lambda3_pred) / lambda3, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(lambda)[3]))
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu1 - mu1_pred) / mu1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(mu)[1]))
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu2 - mu2_pred) / mu2, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(mu)[2]))
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda1 - pred_lambda1_after) / lambda1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="Boost SS      \nError")
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda2 - pred_lambda2_after) / lambda2, color = num_nodes)) + ggplot2::ylim(-1000,5) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda3 - pred_lambda3_after) / lambda3, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu1 - pred_mu1_after) / mu1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu2 - pred_mu2_after) / mu2, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda1 - pred_lambda1_after) / lambda1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="Boost BT      \nError")
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda2 - pred_lambda2_after) / lambda2, color = num_nodes)) + ggplot2::ylim(-1000,5) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda3 - pred_lambda3_after) / lambda3, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu1 - pred_mu1_after) / mu1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu2 - pred_mu2_after) / mu2, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda1 - pred_lambda1_after_lstm_after_dnn) / lambda1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="Boost SS+BT      \nError")
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda2 - pred_lambda2_after_lstm_after_dnn) / lambda2, color = num_nodes)) + ggplot2::ylim(-1000,5) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (lambda3 - pred_lambda3_after_lstm_after_dnn) / lambda3, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu1 - pred_mu1_after_lstm_after_dnn) / mu1, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(durspec, (mu2 - pred_mu2_after_lstm_after_dnn) / mu2, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="MLE Typ      \nError")
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2_r_diff, color = num_nodes)) + ggplot2::ylim(-1000,5) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, mu1_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_typ) + ggplot2::geom_point(ggplot2::aes(durspec, mu2_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, lambda1_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y="MLE Best      \nError")
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, lambda2_r_diff, color = num_nodes)) + ggplot2::ylim(-1000,5) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL, y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, lambda3_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = "Mean duration of speciation",y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, mu1_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)
plot_list_pbd_durspec_rel[[length(plot_list_pbd_durspec_rel) + 1]] <- ggplot2::ggplot(pbd_mle_opt) + ggplot2::geom_point(ggplot2::aes(durspec, mu2_r_diff, color = num_nodes)) + ggplot2::ylim(-10,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale2) + ggplot2::labs(x = NULL,y = NULL)

for (i in 1:length(plot_list_pbd_durspec_rel)) {
  plot_list_pbd_durspec_rel[[i]] <- plot_list_pbd_durspec_rel[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                            legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = x_breaks2, lim = c(0,40)) +
    ggplot2::labs(color = "Tree size")
}

pbd_plot_durspec_rel <- patchwork::wrap_plots(plot_list_pbd_durspec_rel) +
  patchwork::plot_layout(ncol = 5, guides = "collect") +
  patchwork::plot_annotation(title = "Performance Analysis PBD against Mean Duration of Speciation (Relative Error)")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\pbd_plot_durspec_rel.png",
                pbd_plot_durspec_rel,
                device = "png", width = 16,height = 9,dpi = "retina")


## BD Boosting Results
bd_boost_dnn <- readRDS("D:/Habrok/Data/STBO/BD_FREE_TES/BD_FREE_TES_gnn_2_dnn_compensation.rds")
bd_boost_dnn <- bd_boost_dnn %>% dplyr::mutate_all(as.numeric)
bd_boost_dnn <- bd_boost_dnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                          num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                          num_nodes >= 500 ~ "Large"))
bd_boost_dnn_sampled <- bd_boost_dnn[sample(1:nrow(bd_boost_dnn), 2000),]

bd_boost_lstm <- readRDS("D:/Habrok/Data/STBO/BD_FREE_TES/BD_FREE_TES_gnn_2_lstm_compensation.rds")
bd_boost_lstm <- bd_boost_lstm %>% dplyr::mutate_all(as.numeric)
bd_boost_lstm <- bd_boost_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                            num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                            num_nodes >= 500 ~ "Large"))
bd_boost_lstm_sampled <- bd_boost_lstm[sample(1:nrow(bd_boost_lstm), 2000),]

bd_boost_lstm_after_dnn <- readRDS("D:/Habrok/Data/STBO/BD_FREE_TES/BD_FREE_TES_gnn_2_lstm_compensation_after_dnn.rds")
bd_boost_lstm_after_dnn <- bd_boost_lstm_after_dnn %>% dplyr::mutate_all(as.numeric)
bd_boost_lstm_after_dnn <- bd_boost_lstm_after_dnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                                num_nodes >= 500 ~ "Large"))
bd_boost_lstm_after_dnn_sampled <- bd_boost_lstm_after_dnn[sample(1:nrow(bd_boost_lstm_after_dnn), 2000),]

bd_boost_dnn_after_lstm <- readRDS("D:/Habrok/Data/STBO/BD_FREE_TES/BD_FREE_TES_gnn_2_dnn_compensation_after_lstm.rds")
bd_boost_dnn_after_lstm <- bd_boost_dnn_after_lstm %>% dplyr::mutate_all(as.numeric)
bd_boost_dnn_after_lstm <- bd_boost_dnn_after_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                                num_nodes >= 500 ~ "Large"))
bd_boost_dnn_after_lstm_sampled <- bd_boost_dnn_after_lstm[sample(1:nrow(bd_boost_dnn_after_lstm), 2000),]

bd_boost_gnn <- load_final_difference_by_layer(path = "D:/Habrok/Data/STBO",
                                                task_type = "BD_FREE_TES", model = "diffpool", depth = 2)
bd_boost_gnn <- bd_boost_gnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                          num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                          num_nodes >= 500 ~ "Large"))
bd_boost_gnn_sampled <- bd_boost_gnn[sample(1:nrow(bd_boost_gnn), 2000),]

bd_mle_opt <- load_full_mle_result(path = "D:/Habrok/Data/MLE", task_type = "BD_FREE_TES", model_type = "diffpool", no_init = FALSE)

bd_mle_opt <- bd_mle_opt %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                      num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                      num_nodes >= 500 ~ "Large"))

bd_mle_typ <- load_full_mle_result(path = "D:/Habrok/Data/MLE", task_type = "BD_FREE_TES", model_type = "diffpool", no_init = TRUE)

bd_mle_typ <- bd_mle_typ %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                      num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                      num_nodes >= 500 ~ "Large"))

df.borders <- data.frame(intercept=c(0,0,10000),slope=c(1,0,0),Reference=c('zero','same', "mean"))

plot_list_bd_true <- list()
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - lambda_pred, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="GNN      \nError", title = bquote(italic(lambda))) + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - mu_pred, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.72 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL, title = bquote(italic(mu))) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - pred_lambda_after, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="Boost SS      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - pred_mu_after, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.72 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - pred_lambda_after, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="Boost BT      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - pred_mu_after, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.72 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - pred_lambda_after_lstm_after_dnn, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.8 - 0.1) / 2) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="Boost SS+BT      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - pred_mu_after_lstm_after_dnn, color = small)) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (0.72 - 0) / 2)+ ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.72)

plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(lambda, lambda_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="MLE Naive      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(mu, mu_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(lambda, lambda_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y="MLE Best      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(mu, mu_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "True parameter value", y = NULL) + ggplot2::xlim(0,0.72)


for (i in 1:(length(plot_list_bd_true)-4)) {
  plot_list_bd_true[[i]] <- plot_list_bd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list_bd_true)-3):length(plot_list_bd_true)) {
  plot_list_bd_true[[i]] <- plot_list_bd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

bd_plot_true_value <- patchwork::wrap_plots(plot_list_bd_true) +
  patchwork::plot_layout(ncol = 2, guides = "collect", axis_title = "collect_x") + patchwork::plot_annotation(title = "Performance Analysis BD against True Value")

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\bd_plot_true_value3.png",
       bd_plot_true_value,
       width = 10,
       height = 10,
       dpi = "retina",
       device = "png")
ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\bd_plot_true_value3.png",
       bd_plot_true_value,
       width = 10,
       height = 10,
       dpi = "retina",
       device = "png")


bd_mae_gnn  <- bd_boost_gnn_sampled  %>%
  dplyr::group_by(small) %>%
  dplyr::summarise(mae_lambda_gnn = mean(abs(lambda - lambda_pred)),
                   mae_mu_gnn     = mean(abs(mu     - mu_pred)))

bd_mae_boost_dnn <- bd_boost_dnn_sampled %>%
  dplyr::group_by(small) %>%
  dplyr::summarise(mae_lambda_boost_dnn = mean(abs(lambda - pred_lambda_after)),
                   mae_mu_boost_dnn     = mean(abs(mu     - pred_mu_after)))

bd_mae_boost_lstm <- bd_boost_lstm_sampled %>%
  dplyr::group_by(small) %>%
  dplyr::summarise(mae_lambda_boost_lstm = mean(abs(lambda - pred_lambda_after)),
                   mae_mu_boost_lstm     = mean(abs(mu     - pred_mu_after)))

bd_mae_boost_lstm_dnn <- bd_boost_lstm_after_dnn_sampled %>%
  dplyr::group_by(small) %>%
  dplyr::summarise(mae_lambda_boost_lstm_after_dnn = mean(abs(lambda - pred_lambda_after_lstm_after_dnn)),
                   mae_mu_boost_lstm_after_dnn     = mean(abs(mu     - pred_mu_after_lstm_after_dnn)))

bd_mae_mle <- bd_mle_typ %>%
  dplyr::group_by(small) %>%
  dplyr::summarise(mae_lambda_mle = mean(abs(lambda_a_diff)),
                   mae_mu_mle     = mean(abs(mu_a_diff)))

bd_mae_mle_opt <- bd_mle_opt %>%
  dplyr::group_by(small) %>%
  dplyr::summarise(mae_lambda_mle_opt = mean(abs(lambda_a_diff)),
                   mae_mu_mle_opt     = mean(abs(mu_a_diff)))

combined_mae_bd <- bd_mae_gnn %>%          # left‑join the six tables
  dplyr::left_join(bd_mae_boost_dnn,       by = "small") %>%
  dplyr::left_join(bd_mae_boost_lstm,      by = "small") %>%
  dplyr::left_join(bd_mae_boost_lstm_dnn,  by = "small") %>%
  dplyr::left_join(bd_mae_mle,             by = "small") %>%
  dplyr::left_join(bd_mae_mle_opt,         by = "small") %>%
  dplyr::mutate(x_coord = dplyr::case_when(small == "Small"  ~ 100,
                                           small == "Medium" ~ 350,
                                           small == "Large"  ~ 1250),
                y_coord = Inf)             # put the label on the top edge


## BD tree size
nord_palette <-  "lumina"
mae_digits_float <- 2
mae_digits_int <- 0
mae_text_size <- 5
mae_vjust <- 1.6
mae_color <- "#666666"
x_breaks <- c(0, 200, 500, 1000, 1500, 2000)
range_color_scale <- c(0,40)

ggthemr::ggthemr("flat")
plot_list_bd_size <- list()
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y="GNN      \nError", title = bquote(italic(lambda)))
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y = NULL, title = bquote(italic(mu)))
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y="Boost SS      \nError")
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y = NULL)
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y="Boost BT      \nError")
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y = NULL)
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after_lstm_after_dnn)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y="Boost SS+BT      \nError")
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after_lstm_after_dnn)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y = NULL)
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(2*num_nodes, lambda_a_diff)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y="MLE Naive      \nError")
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(2*num_nodes, mu_a_diff)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y = NULL)
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(2*num_nodes, lambda_a_diff)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y="MLE Best      \nError")
plot_list_bd_size[[length(plot_list_bd_size) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(2*num_nodes, mu_a_diff)) + ggplot2::ylim(-1,1) + nord::scale_color_nord(nord_palette, discrete = F, trans="sqrt",limits = range_color_scale) + ggplot2::labs(x = "Number of nodes", y = NULL)

for (i in 1:length(plot_list_bd_size)) {
  plot_list_bd_size[[i]] <- plot_list_bd_size[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 200, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 500, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 2000, color = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = x_breaks, lim = c(0,2100)) +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

label_cols <- c("mae_lambda_gnn",  "mae_mu_gnn",
                "mae_lambda_boost_dnn",  "mae_mu_boost_dnn",
                "mae_lambda_boost_lstm", "mae_mu_boost_lstm",
                "mae_lambda_boost_lstm_after_dnn", "mae_mu_boost_lstm_after_dnn",
                "mae_lambda_mle",  "mae_mu_mle",
                "mae_lambda_mle_opt","mae_mu_mle_opt")

for (i in seq_along(plot_list_bd_size)) {
  plot_list_bd_size[[i]] <-
    plot_list_bd_size[[i]] +
      ggplot2::geom_text(
        data    = combined_mae_bd,
        mapping = ggplot2::aes(x = x_coord,
                               y = y_coord,
                               label = sprintf("%.2f", .data[[ label_cols[i] ]])),
        size     = mae_text_size,
        vjust    = mae_vjust,
        colour   = mae_color,
        fontface = "bold",
        inherit.aes = FALSE)
}

bd_plot_tree_size <- patchwork::wrap_plots(plot_list_bd_size) +
  patchwork::plot_layout(ncol = 2, guides = "collect", axis_titles = "collect_x") +
  patchwork::plot_annotation(title = "Performance Analysis BD against Phylogeny Size")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\bd_plot_tree_size.png",
                bd_plot_tree_size,
                device = "png", width = 10,height = 10,dpi = "retina")



# BD empirical results
condamine_mle <- eveGNN::load_empirical_mle_result("D:\\Habrok\\Data\\EMP")

condamine_mle_bd <- condamine_mle$BD

condamine_gnn_bd <- readRDS("D:\\Habrok\\Data\\EMP\\bd_empirical_gnn_2_lstm_result.rds")

clean_text <- function(text) {
  gsub("\\['(.*?)'\\]", "\\1", text)
}

condamine_gnn_bd$family <- sapply(condamine_gnn_bd$family, clean_text)
condamine_gnn_bd$tree <- sapply(condamine_gnn_bd$tree, clean_text)

condamine_gnn_bd <- condamine_gnn_bd %>% rename(Tree = tree, Family = family)

condamine_data_bd <- left_join(condamine_gnn_bd, condamine_mle_bd, by=c("Family", "Tree"))

condamine_data_bd$num_nodes <- as.numeric(condamine_data_bd$num_nodes)

condamine_data_bd2 <- condamine_data_bd %>% mutate(net_div_nn=pred_lambda_after - pred_mu_after,
                                             net_div_mle=lambda - mu,
                                             sp_nn=pred_lambda_after,
                                             ext_nn=pred_mu_after,
                                             sp_mle=lambda,
                                             ext_mle=mu) %>%
  select(-pred_lambda_after,-pred_lambda_before,-pred_mu_after,-pred_mu_before,-lambda,-mu,-df,-conv,-Model)


ggthemr::ggthemr("flat")

color_palette_emp <- "lumina"
range_color_scale <- range(condamine_data_bd2$num_nodes)

pemp_nd_bd <- ggplot(condamine_data_bd2) + ggpointdensity::geom_pointdensity(aes(x = net_div_nn, y = net_div_mle, color = num_nodes)) +
  lims(x=c(0,3),y=c(0,3)) + geom_abline(slope=1,intercept=0,linewidth=1,linetype="dashed",color="grey") +
  labs(x="Boosting BT", y = "MLE", title = "Net diversification rate (λ - μ)") +
  scale_color_nord(palette = color_palette_emp,discrete=FALSE, limits=range_color_scale, breaks=c(200,500,1000), name="Size") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1/1)

pemp_sp_bd <- ggplot(condamine_data_bd2) + ggpointdensity::geom_pointdensity(aes(x = sp_nn, y = sp_mle, color = num_nodes)) +
  lims(x=c(0,3),y=c(0,3)) + geom_abline(slope=1,intercept=0,linewidth=1,linetype="dashed",color="grey") +
  labs(x="Boosting BT", y = "MLE", title = "Speciation rate (λ)") +
  scale_color_nord(palette = color_palette_emp,discrete=FALSE, limits=range_color_scale, breaks=c(200,500,1000), name="Size") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1/1)

pemp_ext_bd <- ggplot(condamine_data_bd2) + ggpointdensity::geom_pointdensity(aes(x = ext_nn, y = ext_mle, color = num_nodes)) +
  lims(x=c(0,3),y=c(0,3)) + geom_abline(slope=1,intercept=0,linewidth=1,linetype="dashed",color="grey") +
  labs(x="Boosting BT", y = "MLE", title = "Extinction rate (μ)") +
  scale_color_nord(palette = color_palette_emp,discrete=FALSE, limits=range_color_scale, breaks=c(200,500,1000), name="Size") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 1/1)

emp_est_plot_bd <- pemp_nd_bd+pemp_sp_bd+pemp_ext_bd+plot_layout(nrow = 1,guides = "collect") +
  plot_annotation(title="Parameter Estimation on Empirical Trees (Birth-Death Scenario), Neural Network against MLE")

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\bd_empirical.png",
       emp_est_plot_bd,
       width=14,
       height=5,
       dpi="retina",
       device="png")

color_palette_emp  <- "lumina"
range_color_scale  <- range(condamine_data2$num_nodes)

# ── helper: density scatter coloured by size ────────────────────────────
scatter_bd <- function(df, x, y, ttl, lims_xy) {
  ggplot(df) +
    ggpointdensity::geom_pointdensity(
      aes(x = .data[[x]], y = .data[[y]], colour = num_nodes)
    ) +
    lims(x = lims_xy, y = lims_xy) +
    geom_abline(slope = 1, intercept = 0, colour = "grey",
                linewidth = 1, linetype = "dashed") +
    labs(x = "Boosting BT", y = "MLE", title = ttl) +
    scale_color_nord(palette = color_palette_emp, discrete = FALSE,
                     limits  = range_color_scale,
                     breaks  = c(200, 500, 1000),
                     name    = "Size") +
    theme(plot.background  = element_blank(),
          panel.background = element_blank(),
          aspect.ratio     = 1)
}

# ── helper: |Δ|  vs  size panels ───────────────────────────────────────
err_vs_size <- function(df, nn_col, mle_col, ttl, ylim) {
  ggplot(df, aes(x = num_nodes,
                 y = abs(.data[[nn_col]] - .data[[mle_col]]))) +
    geom_point(alpha = .35, size = 1) +
    #geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "black") +
    scale_y_continuous(limits = ylim, name = "|Δ|") +
    scale_x_continuous(name = "Tree size") +
    labs(title = ttl) +
    theme_minimal(base_size = 9) +
    theme(plot.title      = element_text(hjust = .5),
          legend.position = "none",
          aspect.ratio    = .5)
}

# ── top‑row density scatters ────────────────────────────────────────────
pemp_nd_bd  <- scatter_bd(condamine_data2,
                          "net_div_nn", "net_div_mle",
                          "Net diversification rate (λ − μ)", c(0, 3))
pemp_sp_bd  <- scatter_bd(condamine_data2,
                          "sp_nn", "sp_mle",
                          "Speciation rate (λ)", c(0, 3))
pemp_ext_bd <- scatter_bd(condamine_data2,
                          "ext_nn", "ext_mle",
                          "Extinction rate (μ)", c(0, 3))

# ── bottom‑row |Δ| vs size ─────────────────────────────────────────────
e_nd_bd  <- err_vs_size(condamine_data2,
                        "net_div_nn", "net_div_mle",
                        "|λ−μ| vs size", c(0, 3))
e_sp_bd  <- err_vs_size(condamine_data2,
                        "sp_nn", "sp_mle",
                        "|λ| vs size", c(0, 3))
e_ext_bd <- err_vs_size(condamine_data2,
                        "ext_nn", "ext_mle",
                        "|μ| vs size", c(0, 3))

# ── assemble the 2 × 3 layout ──────────────────────────────────────────
emp_est_plot_bd <- (pemp_nd_bd | pemp_sp_bd | pemp_ext_bd) /
  (e_nd_bd   | e_sp_bd   | e_ext_bd) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Parameter Estimation on Empirical Trees (Birth-Death Scenario), Neural Network against MLE"
  )

ggsave("bd_empirical_effect_size.png",
       emp_est_plot_bd,
       width  = 14,
       height = 8,     # taller to accommodate the new second row
       dpi    = "retina",
       device = "png")

## BD parameter distribution
ggplot(bd_boost_dnn) + geom_density(aes(x = lambda))
ggplot(bd_boost_dnn) + geom_density(aes(x = mu))

## Comparison between MLE optimizers, Simplex, Subplex and DEoptim
ddd_simplex_typ <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MLE\\DDD_SIMPLEX", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
ddd_simplex_opt <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MLE\\DDD_SIMPLEX", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

ddd_subplex_typ <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MLE\\DDD_SUBPLEX", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
ddd_subplex_opt <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MLE\\DDD_SUBPLEX", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

ddd_deoptim_typ <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MLE\\DDD_DEOPTIM", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
ddd_deoptim_opt <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MLE\\DDD_DEOPTIM", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

ddd_simplex_typ_new <- load_separated_mle_result(path = "D:\\Habrok\\Data\\STBO-NEW", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
ddd_simplex_opt_new <- load_separated_mle_result(path = "D:\\Habrok\\Data\\STBO-NEW", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

ddd_mmle_typ <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MMLE", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
ddd_mmle_opt <- load_separated_mle_result(path = "D:\\Habrok\\Data\\MMLE", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

ddd_simplex_num_typ <- read_differences(path = "/Users/tianjian/Library/CloudStorage/OneDrive-Persoonlijk/My/Data/NEW", task_type = "DDD", no_init = TRUE)
ddd_simplex_num_opt <- read_differences(path = "/Users/tianjian/Library/CloudStorage/OneDrive-Persoonlijk/My/Data/NEW", task_type = "DDD", no_init = FALSE)
ddd_simplex_num_typ <- convert_differences(ddd_simplex_num_typ,
                    model_type = "DDD_model",
                    task_type   = "DDD_MLE_TES")
ddd_simplex_num_opt <- convert_differences(ddd_simplex_num_opt,
                                           model_type = "DDD_model",
                                           task_type   = "DDD_MLE_TES")
ddd_simplex_ana_typ <- read_differences(path = "/Users/tianjian/Library/CloudStorage/OneDrive-Persoonlijk/My/Data/NEW2", task_type = "DDD", no_init = TRUE)
ddd_simplex_ana_opt <- read_differences(path = "/Users/tianjian/Library/CloudStorage/OneDrive-Persoonlijk/My/Data/NEW2", task_type = "DDD", no_init = FALSE)
ddd_simplex_ana_typ <- convert_differences(ddd_simplex_ana_typ,
                                           model_type = "DDD_model",
                                           task_type   = "DDD_MLE_TES")
ddd_simplex_ana_opt <- convert_differences(ddd_simplex_ana_opt,
                                           model_type = "DDD_model",
                                           task_type   = "DDD_MLE_TES")

filter_and_classify <- function(data, tag) {
  data <-
    data %>% dplyr::mutate(
      small = dplyr::case_when(
        num_nodes < 200 ~ "Small",
        num_nodes >= 200 &
          num_nodes < 500 ~ "Medium",
        num_nodes >= 500 ~ "Large"
      )
    )
  data <-
    data %>% dplyr::filter(
      cap_r_diff < 5000,
      cap_r_diff > -5000,
      lambda_r_diff < 5000,
      lambda_r_diff > -5000,
      mu_r_diff < 5000,
      mu_r_diff > -5000
    )
  data <- data %>% mutate(lambda_pred = lambda - lambda_a_diff, mu_pred = mu - mu_a_diff, cap_pred = cap - cap_a_diff)
  data <- data %>%
    mutate(
      lambda_pred = ifelse(lambda_pred < 0, NA, lambda_pred),
      mu_pred = ifelse(mu_pred < 0, NA, mu_pred),
      cap_pred = ifelse(cap_pred < 0, NA, cap_pred)
    )
  data <- data %>% dplyr::mutate(tag = tag)

  return(data)
}

compute_mae_mle <- function(data) {
  tag = data$tag[1]
  data <- data %>%
    dplyr::group_by(small) %>%
    dplyr::summarize(
      mae_lambda_mle = mean(abs(lambda_a_diff), na.rm = T),
      mae_mu_mle = mean(abs(mu_a_diff), na.rm = T),
      mae_cap_mle = mean(abs(cap_a_diff), na.rm = T),
    ) %>% dplyr::mutate(tag = tag) %>%
  dplyr::mutate(x_coord = dplyr::case_when(small == "Small" ~ 100,
                                           small == "Medium" ~ 350,
                                           small == "Large" ~ 1250)) %>%
  dplyr::mutate(y_coord = Inf)
  return(data)
}

compute_mae_mle2 <- function(data) {
  tag = data$tag[1]
  data <- data %>%
    dplyr::group_by(small) %>%
    dplyr::summarize(
      mae_lambda_mle = mean(abs(lambda-lambda_pred), na.rm = T),
      mae_mu_mle = mean(abs(mu-mu_pred), na.rm = T),
      mae_cap_mle = mean(abs(cap-cap_pred), na.rm = T),
    ) %>% dplyr::mutate(tag = tag) %>%
    dplyr::mutate(x_coord = dplyr::case_when(small == "Small" ~ 100,
                                             small == "Medium" ~ 350,
                                             small == "Large" ~ 1250)) %>%
    dplyr::mutate(y_coord = Inf)
  return(data)
}

ddd_simplex_typ2 <- filter_and_classify(ddd_simplex_typ, "Simplex Typ")
ddd_simplex_opt2 <- filter_and_classify(ddd_simplex_opt, "Simplex Best")

ddd_subplex_typ2 <- filter_and_classify(ddd_subplex_typ, "Subplex Typ")
ddd_subplex_opt2 <- filter_and_classify(ddd_subplex_opt, "Subplex Best")

ddd_deoptim_typ2 <- filter_and_classify(ddd_deoptim_typ, "DEoptim Typ")
ddd_deoptim_opt2 <- filter_and_classify(ddd_deoptim_opt, "DEoptim Best")

ddd_simplex_typ2_new <- filter_and_classify(ddd_simplex_typ_new, "Simplex Typ (NS)")
ddd_simplex_opt2_new <- filter_and_classify(ddd_simplex_opt_new, "Simplex Best (NS)")

ddd_mmle_typ2 <- filter_and_classify(ddd_mmle_typ, "MMLE Typ")
ddd_mmle_opt2 <- filter_and_classify(ddd_mmle_opt, "MMLE Best")

ddd_simplex_num_typ2 <- filter_and_classify(ddd_simplex_num_typ, "Simplex Typ (Num)")
ddd_simplex_num_opt2 <- filter_and_classify(ddd_simplex_num_opt, "Simplex Best (Num)")

ddd_simplex_ana_typ2 <- filter_and_classify(ddd_simplex_ana_typ, "Simplex Typ (Ana)")
ddd_simplex_ana_opt2 <- filter_and_classify(ddd_simplex_ana_opt, "Simplex Best (Ana)")

mle_mae_simp_typ <- compute_mae_mle2(ddd_simplex_typ2)
mle_mae_simp_opt <- compute_mae_mle(ddd_simplex_opt2)

mle_mae_subp_typ <- compute_mae_mle(ddd_subplex_typ2)
mle_mae_subp_opt <- compute_mae_mle(ddd_subplex_opt2)

mle_mae_deop_typ <- compute_mae_mle(ddd_deoptim_typ2)
mle_mae_deop_opt <- compute_mae_mle(ddd_deoptim_opt2)

mle_mae_simp_typ_new <- compute_mae_mle2(ddd_simplex_typ2_new)
mle_mae_simp_opt_new <- compute_mae_mle(ddd_simplex_opt2_new)

mle_mae_mmle_typ <- compute_mae_mle2(ddd_mmle_typ2)
mle_mae_mmle_opt <- compute_mae_mle(ddd_mmle_opt2)

mle_mae_simp_num_typ <- compute_mae_mle2(ddd_simplex_num_typ2)
mle_mae_simp_num_opt <- compute_mae_mle(ddd_simplex_num_opt2)

mle_mae_simp_ana_typ <- compute_mae_mle2(ddd_simplex_ana_typ2)
mle_mae_simp_ana_opt <- compute_mae_mle(ddd_simplex_ana_opt2)

# MLE optimizers against true value
df.borders <- data.frame(intercept=c(0,0),slope=c(1,0),Reference=c('zero','same'))

plot_list_mles <- list()
plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_simplex_typ2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="Simplex Naive      \nError", title = bquote(italic(lambda)))

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_typ2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(mu)))

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_typ2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(K)))

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_simplex_opt2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="Simplex Best      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_opt2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_opt2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_subplex_typ2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_subplex_typ2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="Subplex Naive      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_subplex_typ2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_subplex_typ2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_subplex_typ2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_subplex_typ2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_subplex_opt2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_subplex_opt2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="Subplex Best      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_subplex_opt2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_subplex_opt2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_subplex_opt2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_subplex_opt2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_deoptim_typ2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_deoptim_typ2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="DEoptim Naive      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_deoptim_typ2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_deoptim_typ2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_deoptim_typ2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_deoptim_typ2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_deoptim_opt2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_deoptim_opt2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="Deoptim Best      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_deoptim_opt2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_deoptim_opt2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_deoptim_opt2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_deoptim_opt2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_simplex_typ2_new) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2_new, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="Simplex Typ (NS)      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_typ2_new) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2_new, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_typ2_new) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2_new, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_simplex_opt2_new) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2_new, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="Simplex Best (NS)      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_opt2_new) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2_new, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_opt2_new) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2_new, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_num_typ2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_num_typ2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="New Typ (Numerical)      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_num_typ2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_num_typ2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_num_typ2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_num_typ2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_simplex_num_opt2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_num_opt2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="New Best (Numerical)      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_num_opt2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_num_opt2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::labs(x = "True parameter value", y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_num_opt2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_num_opt2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::labs(color="Tree size") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
    ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
    ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_ana_typ2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_ana_typ2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="New Typ (Analytical)      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_ana_typ2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_ana_typ2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_ana_typ2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_ana_typ2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]] <- ggplot2::ggplot(ddd_simplex_ana_opt2) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_ana_opt2, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="New Best (Analytical)      \nError")

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_ana_opt2) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_ana_opt2, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::labs(x = "True parameter value", y = NULL)

plot_list_mles[[length(plot_list_mles) + 1]]  <- ggplot2::ggplot(ddd_simplex_ana_opt2) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small)) + geom_rug(data = dplyr::filter(ddd_simplex_ana_opt2, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0)))+ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL)


for (i in 1:length(plot_list_mles)) {
  plot_list_mles[[i]] <- plot_list_mles[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

plot_mles_true_value <- patchwork::wrap_plots(plot_list_mles) +
  patchwork::plot_layout(ncol = 3, guides = "collect") + patchwork::plot_annotation(title = "MLE Optimzer Performance against True Value (DDD)")

ggplot2::ggsave("C:\\Source\\eveGNN\\ddd_plot_mles_true_value.png",
                plot_mles_true_value,
                device = "png", width = 14,height = 20,dpi = "retina")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_mles_true_value2.png",
                plot_mles_true_value,
                device = "png", width = 14,height = 12,dpi = "retina")

# ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\ddd_plot_mles_true_value2.png",
#                 plot_mles_true_value,
#                 device = "png", width = 14,height = 12,dpi = "retina")

# MLE optimizers against tree size

nord_palette <-  "aurora"
mae_digits_float <- 2
mae_digits_int <- 0
mae_text_size <- 5
mae_vjust <- 1.6
mae_color <- "#666666"
x_breaks <- c(0, 200, 500, 1000, 1500, 2000)
range_color_scale <- range(c((ddd_simplex_typ2$lambda - ddd_simplex_typ2$mu)/ddd_simplex_typ2$cap),
                           (ddd_simplex_opt2$lambda - ddd_simplex_opt2$mu)/ddd_simplex_opt2$cap,
                           (ddd_subplex_typ2$lambda - ddd_subplex_typ2$mu)/ddd_subplex_typ2$cap,
                           (ddd_subplex_opt2$lambda - ddd_subplex_opt2$mu)/ddd_subplex_opt2$cap,
                            (ddd_deoptim_typ2$lambda - ddd_deoptim_typ2$mu)/ddd_deoptim_typ2$cap,
                            (ddd_deoptim_opt2$lambda - ddd_deoptim_opt2$mu)/ddd_deoptim_opt2$cap)

ggthemr::ggthemr("flat")
plot_list_mles_tree_size <- list()
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]] <- ggplot2::ggplot(ddd_simplex_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y="Simplex Naive      \nError", title = bquote(italic(lambda))) + ggplot2::geom_text(data = mle_mae_simp_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_simplex_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(mu))) + ggplot2::geom_text(data = mle_mae_simp_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_simplex_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_simplex_typ2, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(K))) + ggplot2::geom_text(data = mle_mae_simp_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_simplex_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3)+ nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y="Simplex Best      \nError") + ggplot2::geom_text(data = mle_mae_simp_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_simplex_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_simp_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_simplex_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_simplex_opt2, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale)  + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_simp_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_subplex_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_subplex_typ2, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y="Subplex Naive      \nError") + ggplot2::geom_text(data = mle_mae_subp_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_subplex_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_subplex_typ2, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_subp_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_subplex_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_subplex_typ2, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_subp_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_subplex_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_subplex_opt2, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y="Subplex Best      \nError") + ggplot2::geom_text(data = mle_mae_subp_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_subplex_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_subplex_opt2, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_subp_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_subplex_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_subplex_opt2, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_subp_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_deoptim_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_deoptim_typ2, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y="DEoptim Naive      \nError") + ggplot2::geom_text(data = mle_mae_deop_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_deoptim_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_deoptim_typ2, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_deop_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_deoptim_typ2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_deoptim_typ2, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_deop_typ, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_deoptim_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_deoptim_opt2, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y="DEoptim Best      \nError") + ggplot2::geom_text(data = mle_mae_deop_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_deoptim_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_deoptim_opt2, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(y = NULL,x="Number of nodes") + ggplot2::geom_text(data = mle_mae_deop_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list_mles_tree_size[[length(plot_list_mles_tree_size) + 1]]  <- ggplot2::ggplot(ddd_deoptim_opt2) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) + geom_rug(data = dplyr::filter(ddd_deoptim_opt2, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = mle_mae_deop_opt, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")


for (i in 1:length(plot_list_mles_tree_size)) {
  plot_list_mles_tree_size[[i]] <- plot_list_mles_tree_size[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 200, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 500, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 2000, color = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = x_breaks, lim = c(0,2100)) +
    ggplot2::labs(color = "(λ - μ) / K") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

plot_mles_num_nodes <- patchwork::wrap_plots(plot_list_mles_tree_size) +
  patchwork::plot_layout(ncol = 3, guides = "collect") +
  patchwork::plot_annotation(title = "MLE Optimzer Performance against Phylogeny Size (DDD)")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\ddd_plot_mles_tree_size2.png",
                plot_mles_num_nodes,
                device = "png", width = 14.2,height = 12,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_mles_tree_size2.png",
                plot_mles_num_nodes,
                device = "png", width = 14.2,height = 12,dpi = "retina")

# BD results true value



















load("D:\\Data\\Empirical Trees\\Condamine2019\\FamilyBirdTrees.RData")
is.ultrametric(FamilyBirdTrees$Columbidae$tree)
test_tre <- FamilyBirdTrees$Columbidae$tree
treestats::crown_age(test_tre)
test_tre2 <- eveGNN::rescale_crown_age(test_tre, 10)
treestats::crown_age(test_tre2)
ape::write.nexus(test_tre2, file="D:\\Data\\Empirical Trees\\Condamine2019/trees/tree.tre")
ape::write.nexus(test_tre, file="D:\\Data\\Empirical Trees\\Condamine2019/trees/tree.tre")

test_obj <- parameter_estimation(file_path = "D:/Data/Empirical Trees/Condamine2019/trees",
                     venv_path = "C:/venv/evonn",
                     scenario = "BD")

test_df <- readxl::read_excel("C:\\Users\\tianj\\OneDrive\\桌面\\PlotAP2.xlsx")
test_df_long <- test_df %>% gather(Trait, Value, -site)
ggplot(test_df_long) + geom_boxplot(aes(x=site,y=Value))+facet_wrap(.~Trait, ncol = 2, scales = "free_y")


medium_data <- c(969, 325, 1575, 1216, 566, 8059, 704, 66442, 9122, 317,
                 89, 4539, 1076, 1628, 23528, 1464, 1138, 1390, 302, 223,
                 8777, 193, 1331, 1218, 5449, 17974, 1554, 1061, 853, 293,
                 29873, 3258, 339, 774, 122412, 1264, 601, 377, 816, 1565,
                 3805, 84, 672, 2414, 1880, 731, 607, 135, 586, 551, 642,
                 1199, 2495, 764, 12212, 104, 711, 99, 827, 1316, 205, 425,
                 1439, 2849, 573, 253, 916, 6430, 1003, 634, 824, 705, 1133,
                 404, 724, 1115, 515, 4581, 165, 477, 2187, 366, 404, 16629,
                 1020, 992, 794, 1504, 2957, 170, 116894, 859, 538, 137,
                 437, 5005, 376, 1932, 125, 671)

high_data <-c(649, 6310, 3607, 255, 1072, 37451, 2544, 15711, 1120, 6823, 1085,
              578, 3956, 675, 2769, 73, 2416, 1341, 4428, 1161, 939, 820, 29412,
              1445, 1609, 1066, 1279, 650, 932, 4158, 645, 253, 901, 2276, 2208,
              1217, 500, 287, 3020, 3530, 1109, 268, 271, 2165, 1554, 1326, 9997,
              472, 983, 5320, 8592, 348, 72, 537, 2880, 837, 3199, 1228, 1069,
              496, 5189, 3250, 780, 1509, 1761, 212, 3846, 2118, 251, 163, 1066,
              4268, 375, 1233, 495, 3710, 155, 3964, 952, 524, 1399, 3178, 1360,
              865,  2258, 1773, 1116, 3391, 7617, 732, 495, 9846, 562, 416,
              7557,  1935,  2074,  2093,   203,   705)


example_data <- data.frame(set_num = rep(c("Medium", "High"), each = 100),
                           alllinks = c(medium_data, high_data))

example_data %>% group_by(set_num) %>% summarise(mean(alllinks))

p <- ggplot(example_data, aes(x = set_num, y = alllinks, fill = set_num)) +
  geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "red") +
  coord_trans(y = "log10") +
  scale_fill_brewer(palette = "Blues") +
  theme_classic()
p





# Results Stacking-Boosting
test_ddd <- load_final_difference_by_layer(path = "D:/Habrok/Data/34032024",
                                           task_type = "DDD_FREE_TES", model = "diffpool", depth = 1)

test_ddd <- test_ddd %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                num_nodes >= 500 ~ "Large"))
test_ddd_sampled <- test_ddd[sample(1:nrow(test_ddd), 2000),]

# test_ddd_mle <- load_separated_mle_result(path = "D:/Habrok/Data/MLE/DDD_SIMPLEX", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
# test_ddd_mle_opt <- load_separated_mle_result(path = "D:/Habrok/Data/MLE/DDD_SIMPLEX", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

# Load the new odeint results
test_ddd_mle <- load_separated_mle_result2(path = "C:/Users/tianj/OneDrive/My/Data/COMP3", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
test_ddd_mle_opt <- load_separated_mle_result2(path = "C:/Users/tianj/OneDrive/My/Data/COMP3", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

test_ddd_mle <- test_ddd_mle %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                        num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                        num_nodes >= 500 ~ "Large"))
test_ddd_mle <- test_ddd_mle %>% dplyr::filter(cap_r_diff < 5000, cap_r_diff > -5000,
                                               lambda_r_diff < 5000, lambda_r_diff > -5000,
                                               mu_r_diff < 5000, mu_r_diff > -5000)

test_ddd_mle <- test_ddd_mle %>% mutate(lambda_pred = lambda - lambda_a_diff, mu_pred = mu - mu_a_diff, cap_pred = cap - cap_a_diff)

test_ddd_mle <- test_ddd_mle %>%
  mutate(
    lambda_pred = ifelse(lambda_pred < 0, NA, lambda_pred),
    mu_pred = ifelse(mu_pred < 0, NA, mu_pred),
    cap_pred = ifelse(cap_pred < 0, NA, cap_pred)
  )

test_ddd_mle_opt <- test_ddd_mle_opt %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                num_nodes >= 500 ~ "Large"))
test_ddd_mle_opt <- test_ddd_mle_opt %>% dplyr::filter(cap_r_diff < 5000, cap_r_diff > -5000,
                                                       lambda_r_diff < 5000, lambda_r_diff > -5000,
                                                       mu_r_diff < 5000, mu_r_diff > -5000)

test_ddd_mle_opt <- test_ddd_mle_opt %>% mutate(lambda_pred = lambda - lambda_a_diff, mu_pred = mu - mu_a_diff, cap_pred = cap - cap_a_diff)

test_ddd_mle_opt <- test_ddd_mle_opt %>%
  mutate(
    lambda_pred = ifelse(lambda_pred < 0, NA, lambda_pred),
    mu_pred = ifelse(mu_pred < 0, NA, mu_pred),
    cap_pred = ifelse(cap_pred < 0, NA, cap_pred)
  )

# Bagging results
test_bagging <- readRDS("D:/Habrok/Data/Bagging/DDD_FREE_TES/DDD_FREE_TES_diffpool_1_data.rds")
# convert all columns to numeric
test_bagging <- test_bagging %>% dplyr::mutate_all(as.numeric)
# multiply all columns which name contains "cap" by 1000 using grepl
test_bagging <- test_bagging %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000)

test_bagging<- test_bagging %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                       num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                       num_nodes >= 500 ~ "Large"))

test_bagging_sampled <- test_bagging[sample(1:nrow(test_bagging), 2000),]

test_boost_dnn <- readRDS("D:/Habrok/Data/STBO/DDD_FREE_TES/DDD_FREE_TES_gnn_2_dnn_compensation.rds")
test_boost_dnn <- test_boost_dnn %>% dplyr::mutate_all(as.numeric)
test_boost_dnn <- test_boost_dnn %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000)

test_boost_dnn <- test_boost_dnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                            num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                            num_nodes >= 500 ~ "Large"))
test_boost_dnn_sampled <- test_boost_dnn[sample(1:nrow(test_boost_dnn), 2000),]

test_boost_lstm <- readRDS("D:/Habrok/Data/STBO/DDD_FREE_TES/DDD_FREE_TES_gnn_2_lstm_compensation.rds")
test_boost_lstm <- test_boost_lstm %>% dplyr::mutate_all(as.numeric)
test_boost_lstm <- test_boost_lstm %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000)
test_boost_lstm <- test_boost_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                              num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                              num_nodes >= 500 ~ "Large"))
test_boost_lstm_sampled <- test_boost_lstm[sample(1:nrow(test_boost_lstm), 2000),]

test_boost_lstm_after_dnn <- readRDS("D:/Habrok/Data/STBO/DDD_FREE_TES/DDD_FREE_TES_gnn_2_lstm_compensation_after_dnn.rds")
test_boost_lstm_after_dnn <- test_boost_lstm_after_dnn %>% dplyr::mutate_all(as.numeric)
test_boost_lstm_after_dnn <- test_boost_lstm_after_dnn %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000)
test_boost_lstm_after_dnn <- test_boost_lstm_after_dnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                                  num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                                  num_nodes >= 500 ~ "Large"))
test_boost_lstm_after_dnn_sampled <- test_boost_lstm_after_dnn[sample(1:nrow(test_boost_lstm_after_dnn), 2000),]

test_boost_dnn_after_lstm <- readRDS("D:/Habrok/Data/STBO/DDD_FREE_TES/DDD_FREE_TES_gnn_2_dnn_compensation_after_lstm2.rds")
test_boost_dnn_after_lstm <- test_boost_dnn_after_lstm %>% dplyr::mutate_all(as.numeric)
test_boost_dnn_after_lstm <- test_boost_dnn_after_lstm %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000)
test_boost_dnn_after_lstm <- test_boost_dnn_after_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                                  num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                                  num_nodes >= 500 ~ "Large"))
test_boost_dnn_after_lstm_sampled <- test_boost_dnn_after_lstm[sample(1:nrow(test_boost_dnn_after_lstm), 2000),]

test_boost_gnn <- load_final_difference_by_layer(path = "D:/Habrok/Data/STBO",
                                                 task_type = "DDD_FREE_TES", model = "diffpool", depth = 2)

# Load new results of CNN1D model
test_cnn1d <- load_final_difference_cnn1d(path = "D:/Habrok/Data/STBO",
                                             task_type = "DDD_FREE_TES", model_type = "cnn1d")
# Alternatively, load results from a cloud path
test_cnn1d <- load_final_difference_cnn1d(path = "C:\\Users\\tianj\\OneDrive\\Habrok/Data/STBO",
                                          task_type = "DDD_FREE_TES", model_type = "cnn1d")

test_cnn1d <- test_cnn1d %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                num_nodes >= 500 ~ "Large"))

test_cnn1d_sampled <- test_cnn1d[sample(1:nrow(test_cnn1d), 2000),]

point_alpha <- 0.7
rug_alpha <- 1
df.borders <- data.frame(intercept=c(0,0,10000),slope=c(1,0,0),Reference=c('zero','same', "mean"))

plot_list2 <- list()
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - cap_pred_gnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "GNN")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - cap_pred_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "DNN")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - cap_pred_lstm, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "LSTM")

# plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - cap_pred_median, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "Median")
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_cnn1d_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - cap_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "CNN1D")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - cap_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "Stack")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - pred_cap_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "Boost SS")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - pred_cap_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "Boost BT")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - pred_cap_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "Boost SS+BT")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(cap, cap - pred_cap_after_dnn_after_lstm, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1000 - 10) / 2, show.legend = F)+ ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "Boost BT+SS")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small), alpha = point_alpha) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-1000,1000) + ggplot2::labs(x = NULL, y = NULL, title = "MLE Naive")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(cap, cap-cap_pred, color = small), alpha = point_alpha) +
  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(cap_pred)), aes(cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-1000,1000) + ggplot2::labs(x = "Carrying capactiy K", y=NULL, title = "MLE Best")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot() + 
  ggplot2::geom_tile(aes(x = 500, y = 1000,height=2000,width=Inf), fill="skyblue" , alpha = 0.1) + 
  ggplot2::geom_tile(aes(x = 500, y = -1000,height=2000,width=Inf), fill="salmon", alpha=0.1) + 
  geomtextpath::geom_textabline(slope=1, intercept = -480, label = as.character(expression(hat(y)==bar(y))),
                                lty = "dashed",color = "red", parse=T, hjust = 0.8, size = 5) +
  geomtextpath::geom_textabline(slope=1, intercept = 0, label = as.character(expression(hat(y)==0)),
                                lty = "dotted",color = "purple", parse=T, hjust = 0.8, size = 5) +
  geomtextpath::geom_textabline(slope=0, intercept = 0, label = as.character(expression(hat(y)==y)),
                                lty = "twodash",color = "black", parse=T, hjust = 0.8, size = 5) +
  ggplot2::coord_cartesian(xlim = c(0,1000),ylim=c(-1000,1000))+
  ggplot2::labs(x = NULL, y=NULL, title = "Explanation")+
  ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = 500, ymax = 1000)  +
  ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = -1000, ymax = -500)+ 
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
  
  
for (i in 1:(length(plot_list2)-3)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list2)-2):(length(plot_list2)-1)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

# arrows_df <- data.frame(xs = c(0,1000), xe = c(0,1000), ys=c(333,-333),ye=c(85,-85))
# 
# for (i in 1:length(plot_list2)) {
#   plot_list2[[i]] <- plot_list2[[i]] + 
#     ggplot2::geom_segment(data=arrows_df, aes(x = xs, y = ys, xend = xe, yend = ye),
#                           size=2, arrow = arrow(type = "open", length = unit(0.15, "inches")),
#                           color="salmon",alpha=0.9)
# }

ggthemr::ggthemr("flat")

plot_true_value <- patchwork::wrap_plots(plot_list2) +
  patchwork::plot_layout(ncol = 3, guides = "collect") 

blanklabelplot<-ggplot()+labs(y="Error (True - Estimated)")+theme_classic()+ 
  guides(x = "none", y = "none")

plot_true_value <- blanklabelplot+plot_true_value+ plot_layout(widths=c(1,1000))+ patchwork::plot_annotation(title = "Performance Analysis DDD against True Value of Carrying Capacity")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_cap2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_cap_cnn1d2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")


plot_list2 <- list()
plot_list2[[length(plot_list2) + 1]] <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - lambda_pred_gnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "GNN")

plot_list2[[length(plot_list2) + 1]] <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - lambda_pred_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "DNN")

plot_list2[[length(plot_list2) + 1]] <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - lambda_pred_lstm, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "LSTM")

plot_list2[[length(plot_list2) + 1]] <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - lambda_pred_median, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "Median")

plot_list2[[length(plot_list2) + 1]] <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - lambda_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 *(4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "Stack")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - pred_lambda_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 *(4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "Boost SS")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - pred_lambda_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "Boost BT")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - pred_lambda_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "Boost SS+BT")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - pred_lambda_after_dnn_after_lstm, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2, show.legend = F) + ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "Boost BT+SS")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small), alpha = point_alpha) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y=NULL, title = "MLE Naive")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small), alpha = point_alpha) +
  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-3,3) + ggplot2::labs(x = "Speciation rate λ", y=NULL, title = "MLE Best")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot() + 
  ggplot2::geom_tile(aes(x = 2, y = 3,height=6,width=Inf), fill="skyblue" , alpha = 0.1) + 
  ggplot2::geom_tile(aes(x = 2, y = -3,height=6,width=Inf), fill="salmon", alpha=0.1) + 
  geomtextpath::geom_textabline(slope=1, intercept = -2, label = as.character(expression(hat(y)==bar(y))),
                                lty = "dashed",color = "red", parse=T, hjust = 0.3, size = 5) +
  geomtextpath::geom_textabline(slope=1, intercept = 0, label = as.character(expression(hat(y)==0)),
                                lty = "dotted",color = "purple", parse=T, hjust = 0.3, size = 5) +
  geomtextpath::geom_textabline(slope=0, intercept = 0, label = as.character(expression(hat(y)==y)),
                                lty = "twodash",color = "black", parse=T, hjust = 0.8, size = 5) +
  ggplot2::coord_cartesian(xlim = c(0,4),ylim=c(-3,3))+
  ggplot2::labs(x = NULL, y=NULL, title = "Explanation")+
  ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 3)  +
  ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = -3, ymax = -1.5)+ 
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))

for (i in 1:(length(plot_list2)-3)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list2)-2):(length(plot_list2)-1)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

plot_true_value <- patchwork::wrap_plots(plot_list2) +
  patchwork::plot_layout(nrow = 4, guides = "collect", axis_titles = "collect_y")

blanklabelplot<-ggplot()+labs(y="Error (True - Estimated)")+theme_classic()+ 
  guides(x = "none", y = "none")

plot_true_value <- blanklabelplot+plot_true_value+ plot_layout(widths=c(1,1000)) + patchwork::plot_annotation(title = "Performance Analysis DDD against True Value of Speciation Rate")


ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_lambda2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_lambda_cnn1d2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")

plot_list2 <- list()
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - mu_pred_gnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "GNN")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - mu_pred_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "DNN")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - mu_pred_lstm, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "LSTM")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - mu_pred_median, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "Median")
#plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_cnn1d_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - mu_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "CNN1D")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - mu_pred, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "Stack")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - pred_mu_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "Boost SS")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - pred_mu_after, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "Boost BT")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - pred_mu_after_lstm_after_dnn, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "Boost SS+BT")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(mu, mu - pred_mu_after_dnn_after_lstm, color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (1.5 - 0) / 2, show.legend = F)+ ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y = NULL, title = "Boost BT+SS")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small), alpha = point_alpha) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-2,2) + ggplot2::labs(x = NULL, y=NULL, title = "MLE Naive")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(mu, mu-mu_pred, color = small), alpha = point_alpha) +
  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(mu_pred)), aes(mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-2,2) + ggplot2::labs(x = "Extinction rate μ", y=NULL, title = "MLE Best")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot() + 
  ggplot2::geom_tile(aes(x = 1.5, y = 2,height=4,width=Inf), fill="skyblue" , alpha = 0.1) + 
  ggplot2::geom_tile(aes(x = 1.5, y = -2,height=4,width=Inf), fill="salmon", alpha=0.1) + 
  geomtextpath::geom_textabline(slope=1, intercept = -0.75, label = as.character(expression(hat(y)==bar(y))),
                                lty = "dashed",color = "red", parse=T, hjust = 0.8, size = 4) +
  geomtextpath::geom_textabline(slope=1, intercept = 0, label = as.character(expression(hat(y)==0)),
                                lty = "dotted",color = "purple", parse=T, hjust = 0.8, size = 4) +
  geomtextpath::geom_textabline(slope=0, intercept = 0, label = as.character(expression(hat(y)==y)),
                                lty = "twodash",color = "black", parse=T, hjust = 0.8, size = 4) +
  ggplot2::coord_cartesian(xlim = c(0,1.5),ylim=c(-2,2))+
  ggplot2::labs(x = NULL, y=NULL, title = "Explanation")+
  ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = 1, ymax = 2)  +
  ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = -2, ymax = -1)+ 
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))

for (i in 1:(length(plot_list2)-3)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list2)-2):(length(plot_list2)-1)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

plot_true_value <- patchwork::wrap_plots(plot_list2) +
  patchwork::plot_layout(ncol = 3, guides = "collect") 

blanklabelplot<-ggplot()+labs(y="Error (True - Estimated)")+theme_classic()+ 
  guides(x = "none", y = "none")

plot_true_value <- blanklabelplot+plot_true_value+ plot_layout(widths=c(1,1000))+ patchwork::plot_annotation(title = "Performance Analysis DDD against True Value of Extinction Rate")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_mu2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_mu2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")



plot_list2 <- list()
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_gnn-mu_pred_gnn)/cap_pred_gnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "GNN")
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_dnn-mu_pred_dnn)/cap_pred_dnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "DNN")
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_lstm-mu_pred_lstm)/cap_pred_lstm), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "LSTM")
#plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_median-mu_pred_median)/cap_pred_median), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "Median")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_cnn1d_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred-mu_pred)/cap_pred), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "CNN1D")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred-mu_pred)/cap_pred), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "Stack")
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after-pred_mu_after)/pred_cap_after), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "Boost SS")
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after-pred_mu_after)/pred_cap_after), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "Boost BT")
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after_lstm_after_dnn-pred_mu_after_lstm_after_dnn)/pred_cap_after_lstm_after_dnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "Boost SS+BT")
plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after_dnn_after_lstm-pred_mu_after_dnn_after_lstm)/pred_cap_after_dnn_after_lstm), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -0.001, show.legend = F)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "Boost BT+SS")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes((lambda-mu)/cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred - mu_pred)/(cap_pred)), color = small), alpha = point_alpha) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-0.3,0.3) +
  ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x=NULL,y=NULL,title = "MLE Naive")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) +
  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(lambda_pred)), aes((lambda-mu)/cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred- mu_pred)/(cap_pred)), color = small), alpha = point_alpha) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference), show.legend = F) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-0.3,0.3) +
  ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(x="Carrying capacity effect strength",y=NULL,title = "MLE Best")

plot_list2[[length(plot_list2) + 1]]  <- ggplot2::ggplot() + 
  ggplot2::geom_tile(aes(x = 0.015, y = 0.03,height=0.06,width=Inf), fill="skyblue" , alpha = 0.1) + 
  ggplot2::geom_tile(aes(x = 0.015, y = -0.01,height=0.02,width=Inf), fill="salmon", alpha=0.1) + 
  geomtextpath::geom_textabline(slope=1, intercept = -0.005, label = as.character(expression(hat(y)==bar(y))),
                                lty = "dashed",color = "red", parse=T, hjust = 0.8, size = 4) +
  geomtextpath::geom_textabline(slope=1, intercept = 0, label = as.character(expression(hat(y)==0)),
                                lty = "dotted",color = "purple", parse=T, hjust = 0.8, size = 4) +
  geomtextpath::geom_textabline(slope=0, intercept = 0, label = as.character(expression(hat(y)==y)),
                                lty = "twodash",color = "black", parse=T, hjust = 0.8, size = 4) +
  ggplot2::coord_cartesian(xlim = c(0,0.03),ylim=c(-0.01,0.03))+
  ggplot2::labs(x = NULL, y=NULL, title = "Explanation")+
  ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = 0.015, ymax = 0.03)  +
  ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = -0.01, ymax = -0.005)+ 
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))

for (i in 1:(length(plot_list2)-3)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list2)-2):(length(plot_list2)-1)) {
  plot_list2[[i]] <- plot_list2[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

plot_true_value <- patchwork::wrap_plots(plot_list2) +
  patchwork::plot_layout(ncol = 3, guides = "collect") 

blanklabelplot<-ggplot()+labs(y="Error (True - Estimated)")+theme_classic()+ 
  guides(x = "none", y = "none")

plot_true_value <- blanklabelplot+plot_true_value+ plot_layout(widths=c(1,1000))+ patchwork::plot_annotation(title = "Performance Analysis DDD against True Value of Carrying Capacity Effect Strength")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_cap_effect2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_cap_effect_cnn1d2.png",
                plot_true_value,
                device = "png", width = 14.5,height = 12,dpi = "retina")


ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(lambda, lambda - lambda_pred_gnn, color = small), alpha = point_alpha) +
  ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = -1 * (4.0 - 0.1) / 2) +
  ggplot2::ylim(-3,3) + ggplot2::labs(x = NULL, y="GNN      \nError", title = bquote(italic(lambda)))+
  ggplot2::theme(panel.background = ggplot2::element_blank(),
       legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5),
       panel.grid = ggplot2::element_blank()) +
  ggplot2::labs(color = "Tree size") +
  ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
  ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                             xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0) +
  ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))

ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(lambda, lambda-lambda_pred, color = small), alpha = point_alpha) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes(lambda, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dashed", "same"="dashed", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-3,3)

bagging_mae <- test_bagging %>%
  dplyr::group_by(small) %>%
  dplyr::summarize(mae_lambda_gnn = mean(abs(lambda - lambda_pred_gnn)),
                   mae_mu_gnn = mean(abs(mu - mu_pred_gnn)),
                   mae_cap_gnn = mean(abs(cap - cap_pred_gnn)),
                   mae_lambda_dnn = mean(abs(lambda - lambda_pred_dnn)),
                   mae_mu_dnn = mean(abs(mu - mu_pred_dnn)),
                   mae_cap_dnn = mean(abs(cap - cap_pred_dnn)),
                   mae_lambda_lstm = mean(abs(lambda - lambda_pred_lstm)),
                   mae_mu_lstm = mean(abs(mu - mu_pred_lstm)),
                   mae_cap_lstm = mean(abs(cap - cap_pred_lstm)),
                   mae_lambda_mean = mean(abs(lambda - lambda_pred_mean)),
                   mae_mu_mean = mean(abs(mu - mu_pred_mean)),
                   mae_cap_mean = mean(abs(cap - cap_pred_mean)),
                   mae_lambda_median = mean(abs(lambda - lambda_pred_median)),
                   mae_mu_median = mean(abs(mu - mu_pred_median)),
                   mae_cap_median = mean(abs(cap - cap_pred_median)))
stacking_mae <- test_ddd %>%
  dplyr::group_by(small) %>%
  dplyr::summarize(mae_lambda_stack = mean(abs(lambda - lambda_pred)),
                   mae_mu_stack = mean(abs(mu - mu_pred)),
                   mae_cap_stack = mean(abs(cap - cap_pred))
  )

cnn1d_mae <- test_cnn1d %>%
  dplyr::group_by(small) %>%
  dplyr::summarize(mae_lambda_cnn1d = mean(abs(lambda - lambda_pred)),
                   mae_mu_cnn1d = mean(abs(mu - mu_pred)),
                   mae_cap_cnn1d = mean(abs(cap - cap_pred))
  )

mle_mae <- test_ddd_mle %>%
  filter(lambda_pred != -1,
         mu_pred     != -1,
         cap_pred    != -1) %>%
  group_by(small) %>%
  summarize(
    mae_lambda_mle = mean(abs(lambda - lambda_pred), na.rm = TRUE),
    mae_mu_mle     = mean(abs(mu     - mu_pred),     na.rm = TRUE),
    mae_cap_mle    = mean(abs(cap    - cap_pred),    na.rm = TRUE)
  )

mle_opt_mae <- test_ddd_mle_opt %>%
  filter(lambda_pred != -1,
         mu_pred     != -1,
         cap_pred    != -1) %>%
  group_by(small) %>%
  summarize(
    mae_lambda_mle_opt = mean(abs(lambda_a_diff), na.rm = TRUE),
    mae_mu_mle_opt     = mean(abs(mu_a_diff),     na.rm = TRUE),
    mae_cap_mle_opt    = mean(abs(cap_a_diff),    na.rm = TRUE)
  )
boost_dnn_mae <- test_boost_dnn %>%
  dplyr::group_by(small) %>%
  dplyr::summarize(mae_lambda_boost_dnn = mean(abs(lambda - pred_lambda_after)),
                   mae_mu_boost_dnn = mean(abs(mu - pred_mu_after)),
                   mae_cap_boost_dnn = mean(abs(cap - pred_cap_after))
  )
boost_lstm_mae <- test_boost_lstm %>%
  dplyr::group_by(small) %>%
  dplyr::summarize(mae_lambda_boost_lstm = mean(abs(lambda - pred_lambda_after)),
                   mae_mu_boost_lstm = mean(abs(mu - pred_mu_after)),
                   mae_cap_boost_lstm = mean(abs(cap - pred_cap_after))
  )
boost_lstm_dnn_mae <- test_boost_lstm_after_dnn %>%
  dplyr::group_by(small) %>%
  dplyr::summarize(mae_lambda_boost_lstm_after_dnn = mean(abs(lambda - pred_lambda_after_lstm_after_dnn)),
                   mae_mu_boost_lstm_after_dnn = mean(abs(mu - pred_mu_after_lstm_after_dnn)),
                   mae_cap_boost_lstm_after_dnn = mean(abs(cap - pred_cap_after_lstm_after_dnn))
  )
boost_dnn_lstm_mae <- test_boost_dnn_after_lstm %>%
  dplyr::group_by(small) %>%
  dplyr::summarize(mae_lambda_boost_dnn_after_lstm = mean(abs(lambda - pred_lambda_after_dnn_after_lstm)),
                   mae_mu_boost_dnn_after_lstm = mean(abs(mu - pred_mu_after_dnn_after_lstm)),
                   mae_cap_boost_dnn_after_lstm = mean(abs(cap - pred_cap_after_dnn_after_lstm))
  ) %>%
  #dplyr::mutate(mae_lambda_boost_dnn_after_lstm = dplyr::if_else(mae_lambda_boost_dnn_after_lstm > 10, "++", as.character(mae_lambda_boost_dnn_after_lstm))) %>%
  #dplyr::mutate(mae_mu_boost_dnn_after_lstm = dplyr::if_else(mae_mu_boost_dnn_after_lstm > 10, "++", as.character(mae_mu_boost_dnn_after_lstm))) %>%
  dplyr::mutate(mae_cap_boost_dnn_after_lstm = dplyr::if_else(mae_cap_boost_dnn_after_lstm > 1000, "++", as.character(mae_cap_boost_dnn_after_lstm)))

combined_mae <- dplyr::left_join(bagging_mae, stacking_mae, by = "small") %>%
  dplyr::left_join(mle_mae, by = "small") %>%
  dplyr::left_join(mle_opt_mae, by = "small") %>%
  dplyr::left_join(boost_dnn_mae, by = "small") %>%
  dplyr::left_join(boost_lstm_mae, by = "small") %>%
  dplyr::left_join(boost_lstm_dnn_mae, by = "small") %>%
  dplyr::left_join(boost_dnn_lstm_mae, by = "small") %>%
  dplyr::left_join(cnn1d_mae, by = "small") %>%
  dplyr::mutate(x_coord = dplyr::case_when(small == "Small" ~ 100,
                                           small == "Medium" ~ 350,
                                           small == "Large" ~ 1250)) %>%
  dplyr::mutate(y_coord = Inf)

# Comparison of GNN, DNN and LSTM
test_bagging_sampled <- test_bagging[sample(1:nrow(test_bagging), 2000),]
test_ddd_sampled <- test_ddd[sample(1:nrow(test_ddd), 2000),]
nord_palette <-  "aurora"
mae_digits_float <- 2
mae_digits_int <- 0
mae_text_size <- 4
mae_vjust <- 1
mae_color <- "#666666"
x_breaks <- c(0, 200, 500, 1000, 1500, 2000)
range_color_scale <- range(c((test_bagging_sampled$lambda - test_bagging_sampled$mu)/test_bagging_sampled$cap),
                           (test_ddd_sampled$lambda - test_ddd_sampled$mu)/test_ddd_sampled$cap,
                           (test_ddd_mle_opt$lambda - test_ddd_mle_opt$mu)/test_ddd_mle_opt$cap)

ggthemr::ggthemr("flat")
plot_list <- list()
plot_list[[length(plot_list) + 1]] <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred_gnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y=expression("GNN      \nError"), title = bquote(italic(lambda))) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_gnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred_gnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(mu))) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_gnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred_gnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(K))) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_gnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3)+ nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("DNN      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale)  + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_dnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y=expression("LSTM      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_lstm, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred_median, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y=expression("Median      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_median, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred_median, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_median, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred_median, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_median, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Stack      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_stack, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_stack, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_stack, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost SS\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_boost_dnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost BT\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_boost_lstm, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after_lstm_after_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost SS+BT\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_lstm_after_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after_lstm_after_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_lstm_after_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after_lstm_after_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_boost_lstm_after_dnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after_dnn_after_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost BT+SS\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_dnn_after_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after_dnn_after_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_dnn_after_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after_dnn_after_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = mae_cap_boost_dnn_after_lstm), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")

plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("MLE Naive\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) +   geom_rug(data = dplyr::filter(test_ddd_mle, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(y = NULL,x=NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")

plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) +   geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("MLE Best\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle_opt, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(y = NULL,x="Number of nodes") + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle_opt, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle_opt, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")


for (i in 1:length(plot_list)) {
  plot_list[[i]] <- plot_list[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 200, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 500, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 2000, color = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = x_breaks, lim = c(0,2100)) +
    ggplot2::labs(color = "(λ - μ) / K") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}


plot_num_nodes <- patchwork::wrap_plots(plot_list) +
  patchwork::plot_layout(ncol = 3, guides = "collect") +
  patchwork::plot_annotation(title = "Performance Analysis DDD against Phylogeny Size")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\ddd_plot_tree_size3.png",
                plot_num_nodes,
                device = "png", width = 15,height = 16,dpi = "retina")
ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_tree_size2.png",
                plot_num_nodes,
                device = "png", width = 15,height = 16,dpi = "retina")

ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes(lambda, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("MLE Typ\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")


library(ggh4x)
library(ggnewscale)
library(nord)
library(ggtext)
library(tidyverse)
library(devtools)
library(patchwork)
load_all()

load("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\data and script\\plot.RData")

plot_list <- list()
plot_list[[length(plot_list) + 1]] <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred_gnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y=expression("GNN      \nError"), title = bquote(italic(lambda))) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_gnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred_gnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL, title = bquote(italic(mu))) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_gnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred_gnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y = NULL, title = bquote(italic(K))) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_gnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3)+ nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("DNN      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale)  + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_dnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y=expression("LSTM      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_lstm, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_cnn1d_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y=expression("CNN1D      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_cnn1d, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_cnn1d_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_cnn1d, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_cnn1d_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_cnn1d, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - lambda_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Stack      \nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_stack, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - mu_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_stack, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - cap_pred, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_stack, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost SS\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_boost_dnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost BT\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_boost_lstm, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after_lstm_after_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost SS+BT\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_lstm_after_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after_lstm_after_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_lstm_after_dnn, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after_lstm_after_dnn, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_boost_lstm_after_dnn, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, lambda - pred_lambda_after_dnn_after_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("Boost BT+SS\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_boost_dnn_after_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, mu - pred_mu_after_dnn_after_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_boost_dnn_after_lstm, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_boost_dnn_after_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(num_nodes, cap - pred_cap_after_dnn_after_lstm, color = (lambda - mu) / cap)) + ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = mae_cap_boost_dnn_after_lstm), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")

plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("MLE Naive\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) +   geom_rug(data = dplyr::filter(test_ddd_mle, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) + ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(y = NULL,x=NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")

plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, lambda-lambda_pred, color = (lambda - mu) / cap)) +   geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(lambda_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) + ggplot2::ylim(-3,3) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL,y=expression("MLE Best\nError")) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_lambda_mle_opt, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, mu-mu_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(mu_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-2,2) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(y = NULL,x="Number of nodes") + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_mu_mle_opt, mae_digits_float))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")
plot_list[[length(plot_list) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) + ggplot2::geom_point(ggplot2::aes(num_nodes * 2, cap-cap_pred, color = (lambda - mu) / cap)) +  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(cap_pred)), aes(num_nodes * 2, color = (lambda - mu) / cap), alpha = rug_alpha, size = 2, show.legend = F) +  ggplot2::ylim(-1000,1000) + nord::scale_color_nord(nord_palette, discrete = F, trans = "log", limits = range_color_scale) + ggplot2::labs(x = NULL, y = NULL) + ggplot2::geom_text(data = combined_mae, ggplot2::aes(x = x_coord, y = y_coord, label = paste0(round(mae_cap_mle_opt, mae_digits_int))), size = mae_text_size, vjust = mae_vjust, color = mae_color, fontface = "bold")


for (i in 1:length(plot_list)) {
  plot_list[[i]] <- plot_list[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 200, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 500, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 2000, color = "red", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = x_breaks, lim = c(0,2100)) +
    ggplot2::labs(color = expression(frac(lambda - mu, K))) +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}


plot_num_nodes <- patchwork::wrap_plots(plot_list) +
  patchwork::plot_layout(ncol = 3, guides = "collect") +
  patchwork::plot_annotation(title = "Performance Analysis DDD against Phylogeny Size")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_treesize_cnn1d2.png",
                plot_num_nodes,
                device = "png", width = 15,height = 16,dpi = "retina")

point_alpha <- 0.5
df.borders <- data.frame(intercept=c(0,0,10000),slope=c(1,0,0),Reference=c('zero','same', "mean"))

plot_list3 <- list()
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_gnn-mu_pred_gnn)/cap_pred_gnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "GNN")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_dnn-mu_pred_dnn)/cap_pred_dnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "DNN")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_lstm-mu_pred_lstm)/cap_pred_lstm), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "LSTM")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_ddd_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred-mu_pred)/cap_pred), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "Stack")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after-pred_mu_after)/pred_cap_after), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "Boost SS")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after-pred_mu_after)/pred_cap_after), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "Boost BT")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after_lstm_after_dnn-pred_mu_after_lstm_after_dnn)/pred_cap_after_lstm_after_dnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "Boost SS+BT")

plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_ddd_mle) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes((lambda-mu)/cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred - mu_pred)/(cap_pred)), color = small), alpha = point_alpha) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-0.3,0.3) +
  ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "MLE Typ")

plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) +
  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(lambda_pred)), aes((lambda-mu)/cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred- mu_pred)/(cap_pred)), color = small), alpha = point_alpha) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-0.3,0.3) +
  ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "MLE Best")

for (i in 1:(length(plot_list3)-2)) {
  plot_list3[[i]] <- plot_list3[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size", x = "True carrying capacity effect", y = "Error") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0)
}

for (i in (length(plot_list3)-1):length(plot_list3)) {
  plot_list3[[i]] <- plot_list3[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size", x = "True carrying capacity effect", y = "Error") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0)
}

ggthemr::ggthemr("flat")

plot_true_value <- patchwork::wrap_plots(plot_list3) +
  patchwork::plot_layout(ncol = 3, guides = "collect", axis_titles = "collect") + patchwork::plot_annotation(title = "Performance Analysis DDD Carrying Capacity Effect (λ - μ) / K")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD Start\\Uncategorized\\draft\\Plot\\ddd_plot_true_value_cap2.png",
                plot_true_value,
                device = "png", width = 14,height = 9,dpi = "retina")

plot_list3 <- list()
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_gnn-mu_pred_gnn)/cap_pred_gnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "GNN")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_dnn-mu_pred_dnn)/cap_pred_dnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "DNN")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_bagging_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred_lstm-mu_pred_lstm)/cap_pred_lstm), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "LSTM")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_cnn1d_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred-mu_pred)/cap_pred), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "CNN1D")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after-pred_mu_after)/pred_cap_after), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "Boost SS")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after-pred_mu_after)/pred_cap_after), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "Boost BT")
plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((pred_lambda_after_lstm_after_dnn-pred_mu_after_lstm_after_dnn)/pred_cap_after_lstm_after_dnn), color = small), alpha = point_alpha) + ggplot2::geom_abline(color = "red", lty = "dashed", slope=1, intercept = 0)+ ggplot2::ylim(-0.3,0.3) + ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "Boost SS+BT")

plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_ddd_mle) +
  geom_rug(data = dplyr::filter(test_ddd_mle, is.na(lambda_pred)), aes((lambda-mu)/cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred - mu_pred)/(cap_pred)), color = small), alpha = point_alpha) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-0.3,0.3) +
  ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "MLE Typ")

plot_list3[[length(plot_list3) + 1]]  <- ggplot2::ggplot(test_ddd_mle_opt) +
  geom_rug(data = dplyr::filter(test_ddd_mle_opt, is.na(lambda_pred)), aes((lambda-mu)/cap, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
  ggplot2::geom_point(ggplot2::aes((lambda-mu)/cap, ((lambda-mu)/cap) - ((lambda_pred- mu_pred)/(cap_pred)), color = small), alpha = point_alpha) +
  ggplot2::labs(color="Tree size") +
  ggnewscale::new_scale_color() +
  ggplot2::geom_abline(data = df.borders, aes(lty = Reference, slope=slope, intercept = intercept, color = Reference)) +
  ggplot2::scale_color_manual("Reference", values = c("zero"="purple", "same"="black", "mean"="red"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y)))) +
  ggplot2::scale_linetype_manual("Reference",values = c("zero"="dotted", "same"="twodash", "mean"="dashed"), labels=c('same'=expression(hat(y)==y),'zero'=expression(hat(y)==0),'mean'=expression(hat(y)==bar(y))))+
  ggplot2::ylim(-0.3,0.3) +
  ggplot2::coord_cartesian(xlim=c(0,0.03),ylim=c(-0.01,0.03)) + ggplot2::labs(title = "MLE Best")

for (i in 1:(length(plot_list3)-2)) {
  plot_list3[[i]] <- plot_list3[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size", x = "True carrying capacity effect", y = "Error") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0)
}

for (i in (length(plot_list3)-1):length(plot_list3)) {
  plot_list3[[i]] <- plot_list3[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                      legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size", x = "True carrying capacity effect", y = "Error") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0)
}

ggthemr::ggthemr("flat")

plot_true_value <- patchwork::wrap_plots(plot_list3) +
  patchwork::plot_layout(ncol = 3, guides = "collect", axis_titles = "collect") + patchwork::plot_annotation(title = "Performance Analysis DDD Carrying Capacity Effect (λ - μ) / K")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\ddd_plot_true_value_cap_effect_cnn1d.png",
                plot_true_value,
                device = "png", width = 14,height = 9,dpi = "retina")

remotes::install_github("EvoLandEco/EvoNN")
library(EvoNN)

path <- "c:/test/Columbiformes.tre"

# Estimating parameters under a birth-death scenario:
result1 <- parameter_estimation(file_path = path, scenario = "BD")

# Estimating parameters under a diversity-dependent diversification scenario:
result2 <- parameter_estimation(file_path = path, scenario = "DDD")

bootstrap <- estimation_bootstrap(result2)

ape::read.tree(path) -> tree

treestats::branching_times(tree) -> brts

ddd_ml_result <- DDD::dd_ML(
  brts = brts,
  btorph = 0,
  soc = 2,
  cond = 1,
  ddmodel = 1,
  num_cycles = 1,
  optimmethod = 'simplex'
)

plots <- par(mfrow=c(2,2))
plots <- plot(density(bootstrap$lambda), main = "Speciation rate λ")
plots <- abline(v=result2$pred_lambda, col="red",lty="dashed")
plots <- text(result2$pred_lambda,.5*par('usr')[4],labels=round(result2$pred_lambda,3), col="blue")
plots <- plot(density(bootstrap$mu), main = "Extinction rate μ")
plots <- abline(v=result2$pred_mu, col="red",lty="dashed")
plots <- text(result2$pred_mu,.5*par('usr')[4],labels=round(result2$pred_mu,3), col="blue")
plots <- plot(density(bootstrap$cap), main = "Carrying capacity K")
plots <- abline(v=result2$pred_cap, col="red",lty="dashed")
plots <- text(result2$pred_cap,.5*par('usr')[4],labels=round(result2$pred_cap,3), col="blue")
plots <- plot(density((bootstrap$lambda - bootstrap$mu) / bootstrap$cap), main = "Carrying capacity effect (λ - μ) / K")
plots <- abline(v=((result2$pred_lambda - result2$pred_mu) / result2$pred_cap), col="red",lty="dashed")
plots <- text((result2$pred_lambda - result2$pred_mu) / result2$pred_cap,.5*par('usr')[4],labels=round((result2$pred_lambda - result2$pred_mu) / result2$pred_cap,5), col="blue")



condamine_mle <- load_empirical_mle_result("D:\\Habrok\\Data\\EMP")

condamine_mle_ddd <- condamine_mle$DDD

condamine_gnn_ddd <- readRDS("D:\\Habrok\\Data\\EMP\\empirical_gnn_2_lstm_result.rds")

clean_text <- function(text) {
  gsub("\\['(.*?)'\\]", "\\1", text)
}

condamine_gnn_ddd$family <- sapply(condamine_gnn_ddd$family, clean_text)
condamine_gnn_ddd$tree <- sapply(condamine_gnn_ddd$tree, clean_text)

condamine_gnn_ddd <- condamine_gnn_ddd %>% rename(Tree = tree, Family = family)

condamine_data <- left_join(condamine_gnn_ddd, condamine_mle_ddd, by=c("Family", "Tree"))

condamine_data$num_nodes <- as.numeric(condamine_data$num_nodes)
condamine_data$pred_cap_before <- condamine_data$pred_cap_before * 1000
condamine_data$pred_cap_after <- condamine_data$pred_cap_after * 1000

condamine_list <- condamine_data %>% filter(num_nodes > 300, num_nodes < 1000, cap < 1000)

### Condamine 2019 Ecology Letters Trees Loading
load("D:\\Data\\Empirical Trees\\Condamine2019/FamilyAmphibiaTrees.Rdata")
load("D:\\Data\\Empirical Trees\\Condamine2019/FamilyBirdTrees.Rdata")
load("D:\\Data\\Empirical Trees\\Condamine2019/FamilyCrocoTurtleTrees.Rdata")
load("D:\\Data\\Empirical Trees\\Condamine2019/FamilyMammalTrees.Rdata")
load("D:\\Data\\Empirical Trees\\Condamine2019/FamilySquamateTrees.Rdata")

condamine_tree_list <- list(Amphibia = FamilyAmphibiaTrees,
                            Bird = FamilyBirdTrees,
                            CrocoTurtle = FamilyCrocoTurtleTrees,
                            Mammal = FamilyMammalTrees,
                            Squamate = FamilySquamateTrees)

family_list <- names(condamine_tree_list)

condamine_phy <- list()

for (i in 1:length(family_list)) {
  family_name <- family_list[i]
  tree_list <- names(condamine_tree_list[[family_name]])
  for (j in 1:length(tree_list)) {
    tree_name <- tree_list[j]
    meta <- c("Family" = family_name, "Tree" = tree_name)
    tree <- condamine_tree_list[[family_name]][[tree_name]]$tree
    tree <- eveGNN::rescale_crown_age(tree, 10)

  }
}

# Heatmap BD model
# ────────────────────────────────────────────────────────────────────────
# Helper: make_point_error()
#   Produces a 2‑row patchwork of scatter plots
#   • x‑axis  : true λ
#   • y‑axis  : true μ
#   • point size = |error|            (or magnitude of signed error)
#   • point colour = sign(error)      (blue = under‑estimate, red = over‑estimate)
# ────────────────────────────────────────────────────────────────────────
make_point_error <- function(dat,
                             lambda_col, mu_col,
                             lambda_pred_col = NULL, mu_pred_col = NULL,
                             lambda_diff_col = NULL, mu_diff_col = NULL,
                             mle       = FALSE,
                             abs_err   = TRUE,
                             title_txt = "",
                             max_pt    = 6) {

  # 1. row‑wise error ----------------------------------------------------
  err_fun <- if (abs_err) abs else identity
  have    <- function(x) all(x %in% names(dat))

  if (mle) {
    stopifnot(have(c(lambda_diff_col, mu_diff_col)))
    dat <- dplyr::mutate(dat,
                         lambda_err = err_fun(.data[[lambda_diff_col]]),
                         mu_err     = err_fun(.data[[mu_diff_col]]),
                         lambda_sign = sign(.data[[lambda_diff_col]]),
                         mu_sign     = sign(.data[[mu_diff_col]]))
  } else {
    stopifnot(have(c(lambda_pred_col, mu_pred_col)))
    dat <- dplyr::mutate(dat,
                         lambda_err  = err_fun(.data[[lambda_col]] - .data[[lambda_pred_col]]),
                         mu_err      = err_fun(.data[[mu_col]]     - .data[[mu_pred_col]]),
                         lambda_sign = sign(.data[[lambda_col]] - .data[[lambda_pred_col]]),
                         mu_sign     = sign(.data[[mu_col]]     - .data[[mu_pred_col]]))
  }

  # palette for sign (only used when abs_err = FALSE)
  sign_pal <- scale_colour_manual(values = c("-1" = "#2166ac",
                                             "0"  = "grey60",
                                             "1"  = "#b2182b"),
                                  guide  = "none")

  # 2. plotting ----------------------------------------------------------
  base_aes <- aes(x = .data[[lambda_col]],
                  y = .data[[mu_col]],
                  size = lambda_err)

  pλ <- ggplot(dat, base_aes) +
    geom_point(aes(colour = factor(lambda_sign)), alpha = .6) +
    scale_size_continuous(range = c(0.1, max_pt),
                          name  = if (abs_err)
                            expression("|" * hat(lambda) - lambda * "|")
                          else
                            expression("|" * hat(lambda) - lambda * "| (mag)")) +
  {if (!abs_err) sign_pal} +
    labs(x = expression(True~lambda), y = expression(True~mu),
         title = paste0(title_txt,
                        if (abs_err) ": λ abs. error" else ": λ signed error")) +
    theme_minimal(base_size = 10) +
    theme(plot.title      = element_text(hjust = .5),
          legend.position = "none")

  pμ <- ggplot(dat,
               aes(x = .data[[lambda_col]],
                   y = .data[[mu_col]],
                   size = mu_err)) +
    geom_point(aes(colour = factor(mu_sign)), alpha = .6) +
    scale_size_continuous(range = c(0.1, max_pt),
                          name  = if (abs_err)
                            expression("|" * hat(mu) - mu * "|")
                          else
                            expression("|" * hat(mu) - mu * "| (mag)")) +
  {if (!abs_err) sign_pal} +
    labs(x = expression(True~lambda), y = NULL,
         title = paste0(title_txt,
                        if (abs_err) ": μ abs. error" else ": μ signed error")) +
    theme_minimal(base_size = 10) +
    theme(plot.title      = element_text(hjust = .5),
          legend.position = "none")

  # return a two‑row patchwork
  pλ / pμ
}

# ────────────────────────────────────────────────────────────────────────
# Helper: make_interp_error()
#   Fits a thin‑plate spline to the error field and draws a smooth surface.
#   Arguments are identical to make_point_error(); only the internals differ.
# ────────────────────────────────────────────────────────────────────────
make_interp_error <- function(dat,
                              lambda_col, mu_col,
                              lambda_pred_col = NULL, mu_pred_col = NULL,
                              lambda_diff_col = NULL, mu_diff_col = NULL,
                              mle       = FALSE,
                              abs_err   = TRUE,
                              title_txt = "",
                              grid_res  = 200,      # resolution of the surface
                              ratio_cap = 0.8) {

  # 1. row‑wise error ----------------------------------------------------
  err_fun <- if (abs_err) abs else identity
  have    <- function(x) all(x %in% names(dat))

  if (mle) {
    stopifnot(have(c(lambda_diff_col, mu_diff_col)))
    dat <- dplyr::mutate(dat,
                         lambda_err = err_fun(.data[[lambda_diff_col]]),
                         mu_err     = err_fun(.data[[mu_diff_col]]))
  } else {
    stopifnot(have(c(lambda_pred_col, mu_pred_col)))
    dat <- dplyr::mutate(dat,
                         lambda_err = err_fun(.data[[lambda_col]] - .data[[lambda_pred_col]]),
                         mu_err     = err_fun(.data[[mu_col]]     - .data[[mu_pred_col]]))
  }

  dat <- dplyr::mutate(dat,
                       lambda_s = (lambda - mean(lambda)) / sd(lambda),
                       mu_s     = (mu     - mean(mu))     / sd(mu))

  # 2. helper: interpolate one error column -----------------------------
  interp_surface <- function(zcol) {
    # thin‑plate spline
    fit  <- fields::QTps(as.matrix(dat[, c("lambda_s","mu_s")]), dat[[zcol]], Niterations = 50)

    gx <- seq(min(dat[[lambda_col]]), max(dat[[lambda_col]]), length.out = grid_res)
    gy <- seq(min(dat[[mu_col]]),     max(dat[[mu_col]]),     length.out = grid_res)
    grid <- expand.grid(lambda = gx, mu = gy)

    grid$err <- stats::predict(fit,
                               xnew = as.matrix(grid[, c("lambda_s","mu_s")]))

    # blank biologically impossible area (μ > 0.8 λ)
    grid$err[grid$mu > ratio_cap * grid$lambda] <- NA_real_
    grid
  }

  grid_λ <- interp_surface("lambda_err")
  grid_μ <- interp_surface("mu_err")

  # 3. palette & labels --------------------------------------------------
  if (abs_err) {
    pal_fun  <- function(name)
      scale_fill_nord(palette = "aurora", discrete = FALSE, name = name)
    lab_λ <- expression("|" * hat(lambda) - lambda * "|")
    lab_μ <- expression("|" * hat(mu)     - mu     * "|")
  } else {
    pal_fun <- function(name)
      scale_fill_gradient2(low  = "#2166ac", mid = "white",
                           high = "#b2182b", midpoint = 0, name = name)
    lab_λ <- expression(hat(lambda) - lambda)
    lab_μ <- expression(hat(mu)     - mu)
  }

  # 4. plotting ----------------------------------------------------------
  pλ <- ggplot(grid_λ, aes(lambda, mu, fill = err)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(ambda, mu, z = err), colour = "black", alpha = .4,
                 na.rm = TRUE, inherit.aes = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    pal_fun(lab_λ) +
    labs(x = expression(True~lambda), y = expression(True~mu),
         title = paste0(title_txt,
                        if (abs_err) ": λ abs. error" else ": λ signed error")) +
    coord_equal() +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = .5),
          legend.position = "none")

  pμ <- ggplot(grid_μ, aes(lambda, mu, fill = err)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(ambda, mu, z = err), colour = "black", alpha = .4,
                 na.rm = TRUE, inherit.aes = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    pal_fun(lab_μ) +
    labs(x = expression(True~lambda), y = NULL,
         title = paste0(title_txt,
                        if (abs_err) ": μ abs. error" else ": μ signed error")) +
    coord_equal() +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = .5),
          legend.position = "none")

  # 5. two‑row patchwork -------------------------------------------------
  pλ / pμ
}


plot_list_bd_points_noabs <- list(
  make_interp_error(bd_boost_gnn_sampled,
                   "lambda", "mu",
                   "lambda_pred", "mu_pred",
                   abs_err = FALSE, title_txt = "GNN"),

  make_interp_error(bd_boost_dnn_sampled,
                   "lambda", "mu",
                   "pred_lambda_after", "pred_mu_after",
                   abs_err = FALSE, title_txt = "Boost‑SS"),

  make_interp_error(bd_boost_lstm_sampled,
                   "lambda", "mu",
                   "pred_lambda_after", "pred_mu_after",
                   abs_err = FALSE, title_txt = "Boost‑BT"),

  make_interp_error(bd_boost_lstm_after_dnn_sampled,
                   "lambda", "mu",
                   "pred_lambda_after_lstm_after_dnn",
                   "pred_mu_after_lstm_after_dnn",
                   abs_err = FALSE, title_txt = "Boost‑SS + BT"),

  make_interp_error(bd_mle_typ,
                   "lambda", "mu",
                   lambda_diff_col = "lambda_a_diff",
                   mu_diff_col     = "mu_a_diff",
                   mle = TRUE, abs_err = FALSE, title_txt = "MLE (Typ)"),

  make_interp_error(bd_mle_opt,
                   "lambda", "mu",
                   lambda_diff_col = "lambda_a_diff",
                   mu_diff_col     = "mu_a_diff",
                   mle = TRUE, abs_err = FALSE, title_txt = "MLE (Best)")
)

plot_list_bd_points_abs <- list(
  make_interp_error(bd_boost_gnn_sampled,
                   "lambda", "mu",
                   "lambda_pred", "mu_pred",
                   title_txt = "GNN"),

  make_interp_error(bd_boost_dnn_sampled,
                   "lambda", "mu",
                   "pred_lambda_after", "pred_mu_after",
                   title_txt = "Boost‑SS"),

  make_interp_error(bd_boost_lstm_sampled,
                   "lambda", "mu",
                   "pred_lambda_after", "pred_mu_after",
                   title_txt = "Boost‑BT"),

  make_interp_error(bd_boost_lstm_after_dnn_sampled,
                   "lambda", "mu",
                   "pred_lambda_after_lstm_after_dnn",
                   "pred_mu_after_lstm_after_dnn",
                   title_txt = "Boost‑SS + BT"),

  make_interp_error(bd_mle_typ,
                   "lambda", "mu",
                   lambda_diff_col = "lambda_a_diff",
                   mu_diff_col     = "mu_a_diff",
                   mle = TRUE, title_txt = "MLE (Typ)"),

  make_interp_error(bd_mle_opt,
                   "lambda", "mu",
                   lambda_diff_col = "lambda_a_diff",
                   mu_diff_col     = "mu_a_diff",
                   mle = TRUE, title_txt = "MLE (Best)")
)

bd_pointmap_abs <- patchwork::wrap_plots(plot_list_bd_points_abs) +
  patchwork::plot_layout(ncol = 3, guides = "collect") +
  patchwork::plot_annotation(
    title = "Point‑wise absolute‑error maps in λ–μ space"
  )

bd_pointmap_noabs <- patchwork::wrap_plots(plot_list_bd_points_noabs) +
  patchwork::plot_layout(ncol = 3, guides = "collect") +
  patchwork::plot_annotation(
    title = "Point‑wise signed‑error maps in λ–μ space"
  )

bd_pointmap <- patchwork::wrap_elements(full = bd_pointmap_abs)  |
  patchwork::wrap_elements(full = bd_pointmap_noabs)

ggsave("bd_pointmap.png", bd_pointmap,
       width = 24, height = 16, dpi = 300)



# BD mode, lambda/mu proportion and errors
plot_list_bd_true <- list()
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, lambda - lambda_pred, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y="GNN      \nError", title = bquote(italic(lambda))) + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, mu - mu_pred, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y = NULL, title = bquote(italic(mu))) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, lambda - pred_lambda_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y="Boost SS      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, mu - pred_mu_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, lambda - pred_lambda_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y="Boost BT      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, mu - pred_mu_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, lambda - pred_lambda_after_lstm_after_dnn, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y="Boost SS+BT      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(mu/lambda, mu - pred_mu_after_lstm_after_dnn, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y = NULL) + ggplot2::xlim(0,0.72)

plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(mu/lambda, lambda_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y="MLE Naive      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(mu/lambda, mu_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(mu/lambda, lambda_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y="MLE Best      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(mu/lambda, mu_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Extinction/speciation ratio", y = NULL) + ggplot2::xlim(0,0.72)


for (i in 1:(length(plot_list_bd_true)-4)) {
  plot_list_bd_true[[i]] <- plot_list_bd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list_bd_true)-3):length(plot_list_bd_true)) {
  plot_list_bd_true[[i]] <- plot_list_bd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

bd_plot_true_value <- patchwork::wrap_plots(plot_list_bd_true) +
  patchwork::plot_layout(ncol = 2, guides = "collect", axis_title = "collect_x") + patchwork::plot_annotation(title = "Performance Analysis BD against Extinction/Speciation Ratio")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\bd_plot_true_value_lambda_mu_ratio.png",
                bd_plot_true_value,
                device = "png", width = 10,height = 10,dpi = "retina")

# BD mode, lambda/mu proportion and errors
plot_list_bd_true <- list()
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, lambda - lambda_pred, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y="GNN      \nError", title = bquote(italic(lambda))) + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_gnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, mu - mu_pred, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y = NULL, title = bquote(italic(mu))) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, lambda - pred_lambda_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y="Boost SS      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, mu - pred_mu_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, lambda - pred_lambda_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y="Boost BT      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, mu - pred_mu_after, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, lambda - pred_lambda_after_lstm_after_dnn, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y="Boost SS+BT      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_boost_lstm_after_dnn_sampled) + ggplot2::geom_point(ggplot2::aes(lambda - mu, mu - pred_mu_after_lstm_after_dnn, color = small)) + ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y = NULL) + ggplot2::xlim(0,0.72)

plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(lambda - mu, lambda_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y="MLE Naive      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_typ) + ggplot2::geom_point(ggplot2::aes(lambda - mu, mu_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y = NULL) + ggplot2::xlim(0,0.72)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(lambda - mu, lambda_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y="MLE Best      \nError") + ggplot2::xlim(0,0.8)
plot_list_bd_true[[length(plot_list_bd_true) + 1]] <- ggplot2::ggplot(bd_mle_opt) + ggplot2::geom_point(ggplot2::aes(lambda - mu, mu_a_diff, color = small)) + ggplot2::labs(color="Tree size") +
  ggplot2::ylim(-1,1) + ggplot2::labs(x = "Net diversification rate", y = NULL) + ggplot2::xlim(0,0.72)


for (i in 1:(length(plot_list_bd_true)-4)) {
  plot_list_bd_true[[i]] <- plot_list_bd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

for (i in (length(plot_list_bd_true)-3):length(plot_list_bd_true)) {
  plot_list_bd_true[[i]] <- plot_list_bd_true[[i]] + ggplot2::theme(panel.background = ggplot2::element_blank(),
                                                                    legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(color = "Tree size") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "twodash") +
    ggplot2::annotation_custom(grid::textGrob("Underestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)  +
    ggplot2::annotation_custom(grid::textGrob("Overestimate",gp = grid::gpar(col = "#6C7B8B")),
                               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0) +
    ggplot2::theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5))
}

ggthemr::ggthemr("flat")

bd_plot_true_value <- patchwork::wrap_plots(plot_list_bd_true) +
  patchwork::plot_layout(ncol = 2, guides = "collect", axis_title = "collect_x") + patchwork::plot_annotation(title = "Performance Analysis BD against Net Diversification Rate")

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\bd_plot_true_value_net_div.png",
                bd_plot_true_value,
                device = "png", width = 10,height = 10,dpi = "retina")


## BD Misspecification Results with DDD models
bd_misspec_mle_typ <- load_separated_mle_misspec(path = "D:\\Habrok\\Data\\STBO\\BD_MISSPEC", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
bd_misspec_mle_opt <- load_separated_mle_misspec(path = "D:\\Habrok\\Data\\STBO\\BD_MISSPEC", task_type = "DDD", model_type = "diffpool", no_init = FALSE)
bd_misspec_mle_typ2<- load_separated_mle_misspec(path = "D:\\Habrok\\Data\\STBO\\BD_MISSPEC2", task_type = "DDD", model_type = "diffpool", no_init = TRUE)
bd_misspec_mle_opt2<- load_separated_mle_misspec(path = "D:\\Habrok\\Data\\STBO\\BD_MISSPEC2", task_type = "DDD", model_type = "diffpool", no_init = FALSE)

bd_misspec_mle_typ <- bd_misspec_mle_typ %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                        num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                        num_nodes >= 500 ~ "Large"))
bd_misspec_mle_opt <- bd_misspec_mle_opt %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                        num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                        num_nodes >= 500 ~ "Large"))
bd_misspec_mle_typ2 <- bd_misspec_mle_typ2 %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                          num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                          num_nodes >= 500 ~ "Large"))
bd_misspec_mle_opt2 <- bd_misspec_mle_opt2 %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                            num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                            num_nodes >= 500 ~ "Large"))

bd_misspec_mle_typ <- bd_misspec_mle_typ %>% dplyr::mutate(small = factor(small, levels = c("Large", "Medium", "Small")))
bd_misspec_mle_opt <- bd_misspec_mle_opt %>% dplyr::mutate(small = factor(small, levels = c("Large", "Medium", "Small")))
bd_misspec_mle_typ2 <- bd_misspec_mle_typ2 %>% dplyr::mutate(small = factor(small, levels = c("Large", "Medium", "Small")))
bd_misspec_mle_opt2 <- bd_misspec_mle_opt2 %>% dplyr::mutate(small = factor(small, levels = c("Large", "Medium", "Small")))

bd_misspec_gnn <- load_final_difference_by_layer(path = "D:/Habrok/Data/STBO/BD_MISSPEC",
                               task_type = "DDD_FREE_TES", model = "diffpool", depth = 2)
bd_misspec_gnn <- bd_misspec_gnn %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000)
bd_misspec_gnn <- bd_misspec_gnn %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                                          num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                                          num_nodes >= 500 ~ "Large"))


bd_misspec_gnn_sampled <- bd_misspec_gnn[sample(1:nrow(bd_misspec_gnn), 2000),]

bd_misspec_boost_lstm <- readRDS("D:/Habrok/Data/STBO/BD_MISSPEC/DDD_FREE_TES/DDD_FREE_TES_gnn_2_lstm_compensation.rds")
bd_misspec_boost_lstm <- bd_misspec_boost_lstm %>% dplyr::mutate_all(as.numeric)
bd_misspec_boost_lstm <- bd_misspec_boost_lstm %>% dplyr::mutate_at(dplyr::vars(dplyr::contains("cap")), ~ . * 1000)
bd_misspec_boost_lstm <- bd_misspec_boost_lstm %>% dplyr::mutate(small = dplyr::case_when(num_nodes < 200 ~ "Small",
                                                                              num_nodes >= 200 & num_nodes < 500 ~ "Medium",
                                                                              num_nodes >= 500 ~ "Large"))
bd_misspec_boost_lstm_sampled <- bd_misspec_boost_lstm[sample(1:nrow(bd_misspec_boost_lstm), 2000),]

point_alpha <- 0.4
df.borders <- data.frame(intercept=c(0,0,10000),slope=c(1,0,0),Reference=c('zero','same', "mean"))
theory_df <- tibble(
  div      = seq(min(bd_misspec_mle_typ$lambda - bd_misspec_mle_typ$mu),
                 max(bd_misspec_mle_typ$lambda - bd_misspec_mle_typ$mu),
                 length.out = 200),
  expected = exp(div * 10)
)
props_typ <- bd_misspec_mle_typ %>%
  mutate(R = cap_pred / num_nodes) %>%
  summarise(
    prop_2_5  = mean(R >= 2   & R <  5,        na.rm = TRUE),
    prop_5_20 = mean(R >= 5   & R < 20,        na.rm = TRUE),
    prop_20   = mean(R >= 20  & is.finite(R),  na.rm = TRUE),
    prop_inf  = mean(is.infinite(R))
  )
props_opt <- bd_misspec_mle_opt %>%
  mutate(R = cap_pred / num_nodes) %>%
  summarise(
    prop_2_5  = mean(R >= 2   & R <  5,        na.rm = TRUE),
    prop_5_20 = mean(R >= 5   & R < 20,        na.rm = TRUE),
    prop_20   = mean(R >= 20  & is.finite(R),  na.rm = TRUE),
    prop_inf  = mean(is.infinite(R))
  )
props_typ2 <- bd_misspec_mle_typ2 %>%
    mutate(R = cap_pred / num_nodes) %>%
    summarise(
        prop_2_5  = mean(R >= 2   & R <  5,        na.rm = TRUE),
        prop_5_20 = mean(R >= 5   & R < 20,        na.rm = TRUE),
        prop_20   = mean(R >= 20  & is.finite(R),  na.rm = TRUE),
        prop_inf  = mean(is.infinite(R))
    )
props_opt2 <- bd_misspec_mle_opt2 %>%
    mutate(R = cap_pred / num_nodes) %>%
    summarise(
        prop_2_5  = mean(R >= 2   & R <  5,        na.rm = TRUE),
        prop_5_20 = mean(R >= 5   & R < 20,        na.rm = TRUE),
        prop_20   = mean(R >= 20  & is.finite(R),  na.rm = TRUE),
        prop_inf  = mean(is.infinite(R))
    )
label_text_typ <- with(props_typ, paste0(
  sprintf("2N–5N: %.1f%%", prop_2_5  * 100), "\n",
  sprintf("5N–20N: %.1f%%", prop_5_20 * 100), "\n",
  sprintf("20N–Inf: %.1f%%", prop_20   * 100), "\n",
  sprintf("Inf: %.1f%%", prop_inf  * 100)
))
label_text_opt <- with(props_opt, paste0(
  sprintf("2N–5N: %.1f%%", prop_2_5  * 100), "\n",
  sprintf("5N–20N: %.1f%%", prop_5_20 * 100), "\n",
  sprintf("20N–Inf: %.1f%%", prop_20   * 100), "\n",
  sprintf("Inf: %.1f%%", prop_inf  * 100)
))
label_text_typ2 <- with(props_typ2, paste0(
    sprintf("2N–5N: %.1f%%", prop_2_5  * 100), "\n",
    sprintf("5N–20N: %.1f%%", prop_5_20 * 100), "\n",
    sprintf("20N–Inf: %.1f%%", prop_20   * 100), "\n",
    sprintf("Inf: %.1f%%", prop_inf  * 100)
))
label_text_opt2 <- with(props_opt2, paste0(
    sprintf("2N–5N: %.1f%%", prop_2_5  * 100), "\n",
    sprintf("5N–20N: %.1f%%", prop_5_20 * 100), "\n",
    sprintf("20N–Inf: %.1f%%", prop_20   * 100), "\n",
    sprintf("Inf: %.1f%%", prop_inf  * 100)
))

plot_list_bd_misspec <- list()
plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <- ggplot2::ggplot(bd_misspec_mle_typ) +
  ggplot2::geom_point(ggplot2::aes(lambda - mu, cap_pred, color = small), alpha = point_alpha, inherit.aes = FALSE) +
  ggplot2::guides(colour = "none") +
  ggnewscale::new_scale_color() +
  geom_line(
    data        = theory_df,
    aes(x = div, y = 2 * expected, colour = "Expected nodes"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name   = "Reference",
    values = c("Expected nodes" = "salmon"),
  ) +
  ggplot2::labs(color="Tree size") +
  ggplot2::ylim(0,2000) + ggplot2::labs(x = NULL, y = NULL, title = "MLE Naive") +
  labs(x = NULL, y = "Pseudo-carrying capacity") +
  ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <-
  ggplot(bd_misspec_mle_typ) +
  geom_point(ggplot2::aes(x = num_nodes, y = cap_pred, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,2000) +
    ggplot2::xlim(0,2000) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      name   = "Annotation",
      values = c("K == N" = "gray"),
    ) +
  labs(
    x = NULL,
    y     = NULL
  ) +
    annotate(
      "text",
      x     = Inf, y     = -Inf,
      label = label_text_typ,
      hjust = 1.05, vjust = -0.1,
      color = "black",
      size  = 4
    )+
    ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <- ggplot2::ggplot(bd_misspec_mle_opt) +
  ggplot2::geom_point(ggplot2::aes(lambda - mu, cap_pred, color = small), alpha = point_alpha, inherit.aes = FALSE) +
  ggplot2::guides(colour = "none") +
  ggnewscale::new_scale_color() +
  geom_line(
    data        = theory_df,
    aes(x = div, y = 2 * expected, colour = "Expected nodes"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name   = "Reference",
    values = c("Expected nodes" = "salmon"),
  ) +
  ggplot2::labs(color="Tree size") +
  ggplot2::ylim(0,2000) + ggplot2::labs(x = NULL, y = NULL, title = "MLE Best") +
  labs(x = NULL, y = NULL) +
  ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <-
  ggplot(bd_misspec_mle_opt) +
    geom_point(ggplot2::aes(x = num_nodes, y = cap_pred, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,2000) +
    ggplot2::xlim(0,2000) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      name   = "Annotation",
      values = c("K == N" = "gray"),
    ) +
    labs(
      x = NULL,
      y     = NULL
    ) +
    annotate(
      "text",
      x     = Inf, y     = -Inf,
      label = label_text_opt,
      color = "black",
      hjust = 1.05, vjust = -0.1,
      size  = 4
    ) +
    ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <- ggplot2::ggplot(bd_misspec_mle_typ2) +
  ggplot2::geom_point(ggplot2::aes(lambda - mu, cap_pred, color = small), alpha = point_alpha, inherit.aes = FALSE) +
  ggplot2::guides(colour = "none") +
  ggnewscale::new_scale_color() +
  geom_line(
    data        = theory_df,
    aes(x = div, y = 2 * expected, colour = "Expected nodes"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name   = "Reference",
    values = c("Expected nodes" = "salmon"),
  ) +
  ggplot2::labs(color="Tree size") +
  ggplot2::ylim(0,2000) + ggplot2::labs(x = NULL, y = NULL, title = "MLE Naive") +
  labs(x = NULL, y = "Pseudo-carrying capacity") +
  ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <-
  ggplot(bd_misspec_mle_typ2) +
    geom_point(ggplot2::aes(x = num_nodes, y = cap_pred, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,2000) +
    ggplot2::xlim(0,2000) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      name   = "Annotation",
      values = c("K == N" = "gray"),
    ) +
    labs(
      x = NULL,
      y     = NULL
    ) +
    annotate(
      "text",
      x     = Inf, y     = -Inf,
      label = label_text_typ2,
      hjust = 1.05, vjust = -0.1,
      color = "black",
      size  = 4
    )+
    ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <- ggplot2::ggplot(bd_misspec_mle_opt2) +
  ggplot2::geom_point(ggplot2::aes(lambda - mu, cap_pred, color = small), alpha = point_alpha, inherit.aes = FALSE) +
  ggplot2::guides(colour = "none") +
  ggnewscale::new_scale_color() +
  geom_line(
    data        = theory_df,
    aes(x = div, y = 2 * expected, colour = "Expected nodes"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name   = "Reference",
    values = c("Expected nodes" = "salmon"),
  ) +
  ggplot2::labs(color="Tree size") +
  ggplot2::ylim(0,2000) + ggplot2::labs(x = NULL, y = NULL, title = "MLE Best") +
  labs(x = NULL, y = NULL) +
  ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <-
  ggplot(bd_misspec_mle_opt2) +
    geom_point(ggplot2::aes(x = num_nodes, y = cap_pred, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,2000) +
    ggplot2::xlim(0,2000) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      name   = "Annotation",
      values = c("K == N" = "gray"),
    ) +
    labs(
      x = NULL,
      y     = NULL
    ) +
    annotate(
      "text",
      x     = Inf, y     = -Inf,
      label = label_text_opt2,
      color = "black",
      hjust = 1.05, vjust = -0.1,
      size  = 4
    ) +
    ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <- ggplot2::ggplot(bd_misspec_gnn_sampled) +
  ggplot2::geom_point(ggplot2::aes(lambda - mu, cap_diff, color = small), alpha = point_alpha, inherit.aes = FALSE) +
  ggplot2::guides(colour = "none") +
  ggnewscale::new_scale_color() +
  geom_line(
    data        = theory_df,
    aes(x = div, y = 2 * expected, colour = "Expected nodes"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  geom_hline(
    aes(yintercept = 445, colour = "Conditional mean"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name   = "Reference",
    values = c("Expected nodes" = "salmon", "Conditional mean" = "steelblue")
  ) +
  ggplot2::labs(color="Tree size") +
  ggplot2::ylim(0,2000) + ggplot2::labs(x = NULL, y = NULL, title = "GNN") +
  labs(x = "Net diversification rate (λ - μ)", y = "Pseudo-carrying capacity")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <-
  ggplot(bd_misspec_gnn_sampled) +
    geom_point(ggplot2::aes(x = num_nodes, y = cap_diff, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,2000) +
    ggplot2::xlim(0,2000) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      name   = "Annotation",
      values = c("K == N" = "gray"),
    ) +
    labs(
      x = "Number of nodes",
      y     = NULL
    ) +
    ggplot2::guides(colour = "none")

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <- ggplot2::ggplot(bd_misspec_boost_lstm_sampled) +
  ggplot2::geom_point(ggplot2::aes(lambda - mu, pred_cap_after, color = small), alpha = point_alpha, inherit.aes = FALSE) +
  guides(
    colour = guide_legend(
      override.aes = list(
        alpha   = c(1, 1, 1)
      )
    ),
    title = "Tree size"
  ) +
  ggnewscale::new_scale_color() +
  geom_line(
    data        = theory_df,
    aes(x = div, y = 2 * expected, colour = "Expected nodes"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  geom_hline(
    aes(yintercept = 445, colour = "Conditional mean"),
    linetype    = "dotdash",
    size        = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    name   = "Reference",
    values = c("Expected nodes" = "salmon", "Conditional mean" = "steelblue")
  ) +
  ggplot2::labs(color="Tree size") +
  ggplot2::ylim(0,2000) + ggplot2::labs(x = NULL, y = NULL, title = "Boost BT") +
  labs(x = "Net diversification rate (λ - μ)", y = NULL)

plot_list_bd_misspec[[length(plot_list_bd_misspec) + 1]] <-
  ggplot(bd_misspec_boost_lstm_sampled) +
    geom_point(ggplot2::aes(x = num_nodes, y = pred_cap_after, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,2000) +
    ggplot2::xlim(0,2000) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      name   = "Annotation",
      values = c("K == N" = "gray"),
    ) +
    labs(
      x = "Number of nodes",
      y     = NULL
    ) +
    ggplot2::guides(colour = "none")

bd_plot_misspec <- patchwork::wrap_plots(plot_list_bd_misspec) +
  patchwork::plot_layout(ncol = 4, guides = "collect", axis_title = "collect_x") +
  patchwork::plot_annotation(title = "Estimate BD Trees under DDD Assumption (Carrying Capacity)") &
  ggplot2::theme(
    aspect.ratio =  1/ 1,
    panel.background = ggplot2::element_blank()
  )

ggplot2::ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\bd_plot_misspec_lambda_mu2.png",
                bd_plot_misspec,
                device = "png", width = 15,height = 12,dpi = "retina")

rug_alpha <- 0.4
plot_list_bd_misspec_nd <- list()
plot_list_bd_misspec_nd[[length(plot_list_bd_misspec_nd) + 1]] <-
  ggplot(bd_misspec_mle_typ) +
    geom_point(ggplot2::aes(x = lambda - mu, y = lambda_pred - mu_pred, color = small), alpha = 0.6, size = 2) +
    geom_rug(data = dplyr::filter(bd_misspec_mle_typ, is.na(lambda_pred - mu_pred)), aes(lambda - mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,8) +
    ggplot2::xlim(0,1) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(x = NULL, y = "Pseudo net diversification rate", title = "MLE Typ")


plot_list_bd_misspec_nd[[length(plot_list_bd_misspec_nd) + 1]] <-
  ggplot(bd_misspec_mle_opt) +
    geom_point(ggplot2::aes(x = lambda - mu, y = lambda_pred - mu_pred, color = small), alpha = 0.6, size = 2) +
    geom_rug(data = dplyr::filter(bd_misspec_mle_opt, is.na(lambda_pred - mu_pred)), aes(lambda - mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,8) +
    ggplot2::xlim(0,1) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "MLE Best")

plot_list_bd_misspec_nd[[length(plot_list_bd_misspec_nd) + 1]] <-
  ggplot(bd_misspec_mle_typ2) +
    geom_point(ggplot2::aes(x = lambda - mu, y = lambda_pred - mu_pred, color = small), alpha = 0.6, size = 2) +
    geom_rug(data = dplyr::filter(bd_misspec_mle_typ2, is.na(lambda_pred - mu_pred)), aes(lambda - mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,8) +
    ggplot2::xlim(0,1) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(x = NULL, y = "Pseudo net diversification rate", title = "MLE Typ (Numerical Solver)")

plot_list_bd_misspec_nd[[length(plot_list_bd_misspec_nd) + 1]] <-
  ggplot(bd_misspec_mle_opt2) +
    geom_point(ggplot2::aes(x = lambda - mu, y = lambda_pred - mu_pred, color = small), alpha = 0.6, size = 2) +
    geom_rug(data = dplyr::filter(bd_misspec_mle_opt2, is.na(lambda_pred - mu_pred)), aes(lambda - mu, color = small), alpha = rug_alpha, size = 2, show.legend = F) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,8) +
    ggplot2::xlim(0,1) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "MLE Best (Numerical Solver)")

plot_list_bd_misspec_nd[[length(plot_list_bd_misspec_nd) + 1]] <-
  ggplot(bd_misspec_gnn_sampled) +
    geom_point(ggplot2::aes(x = lambda - mu, y = lambda_pred - mu_pred, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,8) +
    ggplot2::xlim(0,1) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::guides(colour = "none") +
    ggplot2::labs(x = "True net diversification rate (λ - μ)", y = "Pseudo net diversification rate", title = "GNN")

plot_list_bd_misspec_nd[[length(plot_list_bd_misspec_nd) + 1]] <-
  ggplot(bd_misspec_boost_lstm_sampled) +
    geom_point(ggplot2::aes(x = lambda - mu, y = pred_lambda_after - pred_mu_after, color = small), alpha = 0.6, size = 2) +
    ggplot2::guides(colour = "none") +
    ggplot2::ylim(0,8) +
    ggplot2::xlim(0,1) +
    ggnewscale::new_scale_color() +
    geom_abline(aes(color = "K == N", slope = 1, intercept = 0),
                linetype = "dashed"
    ) +
    ggplot2::guides(colour = "none")  +
    ggplot2::labs(x = "True net diversification rate (λ - μ)", y = NULL, title = "Boost BT")

plot_list_bd_misspec_nd <- patchwork::wrap_plots(plot_list_bd_misspec_nd) +
  patchwork::plot_layout(ncol = 2, guides = "collect", axis_title = "collect_x") +
  patchwork::plot_annotation(title = "Estimate BD Trees under DDD Assumption (Net Diversification Rate)") &
  ggplot2::theme(
    aspect.ratio =  1/ 1,
    panel.background = ggplot2::element_blank()
  )

ggplot2::ggsave("bd_plot_misspec_lambda_mu_nd.png",
                plot_list_bd_misspec_nd,
                device = "png", width = 8,height = 12,dpi = "retina")

summarize_ddml_errors <- function(
  err_dir = "safe_ddml_errors",
  markdown = FALSE,
  md_file = "safe_ddml_errors_report.md"
) {
  if (!dir.exists(err_dir)) stop("Directory ‘", err_dir, "’ not found.")
  f <- list.files(err_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(f)) { message("No error files."); return(invisible(NULL)) }

  rows <- lapply(f, \(z) {
    y <- readRDS(z)
    data.frame(
      opt_method = y$opt_method,
      int_method = y$int_method,
      error_type = y$type,
      error_msg  = y$error_msg,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)

  by_opt  <- as.data.frame(table(df$opt_method),  stringsAsFactors = FALSE)
  by_int  <- as.data.frame(table(df$int_method),  stringsAsFactors = FALSE)
  by_type <- as.data.frame(table(df$error_type), stringsAsFactors = FALSE)
  names(by_opt)  <- c("opt_method", "n")
  names(by_int)  <- c("int_method", "n")
  names(by_type) <- c("error_type", "n")

  combo_tally <- as.data.frame(table(df$opt_method, df$int_method, df$error_type),
                               stringsAsFactors = FALSE)
  names(combo_tally) <- c("opt_method", "int_method", "error_type", "n")
  combo_tally <- combo_tally[combo_tally$n > 0, ]

  combo_msgs <- aggregate(
    df$error_msg,
    list(opt_method = df$opt_method,
         int_method = df$int_method,
         error_msg  = df$error_msg),
    length
  )
  names(combo_msgs)[4] <- "count"

  if (markdown) {
    md <- c(
      "# safe_ddml error summary",
      "",
      "## Counts by optimiser",
      knitr::kable(by_opt, format = "markdown"),
      "",
      "## Counts by solver",
      knitr::kable(by_int, format = "markdown"),
      "",
      "## Counts by error type",
      knitr::kable(by_type, format = "markdown"),
      "",
      "## Counts by optimiser × solver × error type",
      knitr::kable(combo_tally, format = "markdown"),
      ""
    )
    split_msgs <- split(combo_msgs, list(combo_msgs$opt_method,
                                         combo_msgs$int_method),
                        drop = TRUE)
    for (nm in names(split_msgs)) {
      part <- split_msgs[[nm]]
      md <- c(
        md,
        sprintf("### %s", gsub("\\.", " / ", nm)),
        knitr::kable(part[, c("error_msg", "count")], format = "markdown"),
        ""
      )
    }
    writeLines(md, md_file)
  }

  list(
    counts_by_optimizer = by_opt,
    counts_by_solver    = by_int,
    counts_by_type      = by_type,
    counts_by_combo     = combo_tally,
    messages_by_combo   = combo_msgs
  )
}

error_path <- "/Users/tianjian/Library/CloudStorage/OneDrive-Persoonlijk/My/Data/DDD_MLE_TES/safe_ddml_errors"
report_path <- "ddml_errors_report.md"

summarize_ddml_errors(error_path, markdown = T, report_path) -> ddml_errors
ddml_errors

get_failed_runs <- function(
  err_dir   = "safe_ddml_errors",
  opt_pick  = NULL,     # e.g. "simplex" or c("simplex","DEoptim") or "all"
  ode_pick  = NULL,     # e.g. "analytical" or "odeint::runge_kutta_cash_karp54" or "all"
  msg_pick  = NULL,     # error-message substring(s) or "all"
  ask_when_missing = interactive()
) {
  if (!dir.exists(err_dir)) stop("Directory ‘", err_dir, "’ not found.")
  rds_files <- list.files(err_dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(rds_files)) stop("No *.rds files in ‘", err_dir, "’.")

  data_tbl <- do.call(
    rbind,
    lapply(rds_files, function(f) {
      x <- readRDS(f)
      data.frame(
        opt_method = x$opt_method,
        ode_solver = x$int_method,
        error_msg  = x$error_msg,
        brts       = I(list(x$brts)),
        initpars   = I(list(x$initpars)),
        stringsAsFactors = FALSE
      )
    })
  )

  ## menu helper --------------------------------------------------
  choose_one <- function(choices, label) {
    pick <- utils::menu(c(choices, "ALL"), title = paste("Select", label))
    if (pick == 0) stop("Selection cancelled.")
    if (pick == length(choices) + 1) "all" else choices[pick]
  }

  ## optimiser selection -----------------------------------------
  if (is.null(opt_pick) && ask_when_missing)
    opt_pick <- choose_one(unique(data_tbl$opt_method), "optimizer")
  if (is.null(opt_pick) || identical(opt_pick, "all"))
    opt_pick <- unique(data_tbl$opt_method)

  data_opt <- subset(data_tbl, opt_method %in% opt_pick)
  if (!nrow(data_opt)) {
    warning("No failures for the chosen optimiser(s).")
    return(list())
  }

  ## solver selection --------------------------------------------
  if (is.null(ode_pick) && ask_when_missing)
    ode_pick <- choose_one(unique(data_opt$ode_solver), "ODE solver")
  if (is.null(ode_pick) || identical(ode_pick, "all"))
    ode_pick <- unique(data_opt$ode_solver)

  data_ode <- subset(data_opt, ode_solver %in% ode_pick)
  if (!nrow(data_ode)) {
    warning("No failures for the chosen optimiser / solver combination.")
    return(list())
  }

  ## message selection -------------------------------------------
  avail_msgs <- unique(data_ode$error_msg)
  if (is.null(msg_pick) && ask_when_missing)
    msg_pick <- choose_one(avail_msgs, "error message")
  if (is.null(msg_pick) || identical(msg_pick, "all"))
    msg_pick <- avail_msgs

  data_final <- subset(data_ode, error_msg %in% msg_pick)
  if (!nrow(data_final)) {
    warning("No failures match the specified error message(s).")
    return(list())
  }

  split(
    data_final[, c("brts", "initpars", "opt_method", "ode_solver", "error_msg")],
    seq_len(nrow(data_final))
  )
}

# Interactive filters (let the menu drive choices)
sel <- get_failed_runs(error_path)

# Or manually retrieve only simplex + analytical failures, all errors
sx <- get_failed_runs(error_path, opt_pick = "simplex", ode_pick = "analytical",
                      msg_pick = "all", ask_when_missing = FALSE)

# ────────────────────────────────────────────────────────────────
# 1.  Helper  ─ read saved differences*.rds into a tidy data frame
# ────────────────────────────────────────────────────────────────
read_differences <- function(path,
                             task_type = c("DDD", "PBD"),
                             no_init   = FALSE) {
  task_type <- match.arg(task_type)

  subdir <- if (no_init) file.path(path, paste0(task_type, "_MLE_TES"), "NO_INIT")
  else         file.path(path, paste0(task_type, "_MLE_TES"))

  files  <- list.files(subdir, pattern = "^differences_[0-9]+\\.rds$", full.names = TRUE)
  if (length(files) == 0)
    stop("No differences_*.rds files found in ", subdir)

  purrr::map_dfr(files, \(f) {
    raw <- readRDS(f)

    # Handle legacy one-run structure vs. new 6-run list
    runs <- if (is.list(raw) && !is.null(raw$opt_method)) list(raw) else raw

    purrr::map_dfr(runs, \(x) {
      if (is.null(x$differences) || is.na(x$loglik)) return(NULL)

      tibble::tibble(
        file        = basename(f),
        opt_method  = x$opt_method,
        integrator  = x$methode,
        loglik      = x$loglik,
        param       = c("lambda", "mu", "K"),
        true        = x$differences$true,
        pred        = x$differences$mle,
        abs_err     = x$differences$true - x$differences$mle,
        rel_err_pc  = 100 * (x$differences$true - x$differences$mle) /
          x$differences$true,
        nnode       = x$nnode
      )
    })
  })
}
# ────────────────────────────────────────────────────────────────
# 2.  Visualisation  ─ log-likelihood vs. estimation error
# ────────────────────────────────────────────────────────────────
plot_ll_vs_error <- function(df,
                             tag = NULL,
                             rel_cut = 2000) {
  library(ggplot2)
  library(ggtext)          # markdown strip labels

  # ── 1 ‒ keep only rows within ±rel_cut % ───────────────────────
  df <- dplyr::filter(df, abs(rel_err_pc) <= rel_cut)

  # reorder facetting vars (optional, for nicer layout)
  df$param      <- factor(df$param, levels = c("lambda", "mu", "K"))
  df$opt_method <- factor(df$opt_method, levels = c("simplex", "subplex", "DEoptim"))

  # ── 2 ‒ plot ───────────────────────────────────────────────────
  p <- ggplot(df,
              aes(x = loglik,
                  y = rel_err_pc,
                  colour = integrator)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_grid(rows = vars(param),
               cols = vars(opt_method),
               scales = "free_y",
               labeller = labeller(param = greek_param_labeller)) +
    scale_colour_brewer(palette = "Dark2",
                        name = "Integrator") +
    labs(
      x = "Maximum log-likelihood  (higher = better fit)",
      y = "Relative error  (%)",
      title = paste0("Estimation Error (≤ ±2000 %) vs. Maximized Log-Likelihood", " (", tag, ")"),
      subtitle = glue::glue("Facetted by Parameter (Rows) and Optimizer (Columns)")
    ) +
    theme(
      legend.position = "bottom",
      plot.background = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      legend.background = element_blank()
    )

  p
}

result_path <- "/Users/tianjian/Library/CloudStorage/OneDrive-Persoonlijk/My/Data/"

mle_accu_typ <- read_differences(result_path, task_type = "DDD", no_init = TRUE)
mle_accu_opt <- read_differences(result_path, task_type = "DDD", no_init = FALSE)

# plot absolute error
p_mle_accu_typ <- plot_ll_vs_error(mle_accu_typ, "Typical Case")
p_mle_accu_opt <- plot_ll_vs_error(mle_accu_opt, "Best Case")

p_mle_accu <- p_mle_accu_typ / p_mle_accu_opt

ggplot2::ggsave("mle_accuracy.png", p_mle_accu, device = "png",
                width = 10, height = 12, dpi = "retina")

library(dplyr)

add_best_of <- function(df) {
  df %>%
    group_by(file, integrator) %>%                    # one data set per integrator
    mutate(best_opt = opt_method[ which.max(loglik) ]) %>%   # optimiser that won
    mutate(strategy = ifelse(opt_method == best_opt,
                             "MMLE", opt_method)) %>%   # label every row
    ungroup()
}


plot_best_vs_single <- function(df,
                                err   = "rel_err_pc",   # or "abs_err"
                                ylim  = c(-100, 100)) { # zoom to ±100 %

  df <- add_best_of(df) %>%
    filter(abs(rel_err_pc) <= 2000)   # ← same 2000 % guard

  # order legend nicely
  df$strategy <- factor(df$strategy,
                        levels = c("simplex", "subplex", "DEoptim", "best-of-3"))

  library(ggplot2)
  ggplot(df, aes(x = strategy,
                 y = .data[[err]],
                 fill = strategy)) +
    geom_violin(trim = TRUE, width = 0.85, alpha = 0.6,
                colour = "grey30", linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3) +
    facet_grid(rows = vars(param),
               cols = vars(integrator),
               scales = "free_y",
               labeller = labeller(param = greek_param_labeller)) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    coord_cartesian(ylim = ylim) +
    labs(
      x = NULL,
      y = if (err == "abs_err") "Absolute error (true – est.)"
      else "Relative error  (%)",
      title = "Does choosing the highest log-likelihood help?",
      subtitle = "“best-of-3” = run with the largest LL among simplex / subplex / DEoptim"
    ) +
    theme(
      strip.text  = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

plot_best_vs_single(df, err = "rel_err_pc", ylim = c(-200, 200))

greek_param_labeller <- ggplot2::as_labeller(
  c(lambda = "lambda",     # parsed to λ
    mu     = "mu",         # parsed to μ
    K      = "italic(K)"), # or just "K"
  default = ggplot2::label_parsed
)

plot_typ_vs_opt <- function(df_typ,
                            df_opt,
                            err      = "rel_err_pc",
                            rel_cut  = 2000,
                            ylim     = c(-100, 100),
                            smooth   = FALSE) {

  df_typ <- add_best_of(df_typ) %>% mutate(case = "Typical")
  df_opt <- add_best_of(df_opt) %>% mutate(case = "Best")

  df <- dplyr::bind_rows(df_typ, df_opt) %>%
    dplyr::filter(abs(rel_err_pc) <= rel_cut)

  df$param    <- factor(df$param, levels = c("lambda", "mu", "K"))
  df$strategy <- factor(df$strategy,
                        levels = c("simplex", "subplex", "DEoptim", "MMLE"))
  df$case     <- factor(df$case, levels = c("Best", "Typical"))

  library(ggplot2)

  pos <- position_dodge(width = 0.75)

  ggplot(df,
         aes(x = strategy,
             y = .data[[err]],
             fill = case)) +
    geom_violin(position = pos, alpha = 0.65, width = 0.9, colour = "grey30") +
    geom_boxplot(position = pos, width = 0.15,
                 outlier.shape = NA, linewidth = 0.3) +
  {if (smooth)
    geom_smooth(aes(colour = case, group = interaction(case, strategy)),
                method = "loess", se = FALSE, linewidth = 0.5,
                position = pos)} +
    facet_grid(rows = vars(param), cols = vars(integrator), scales = "free_y",
               labeller = labeller(param = greek_param_labeller)) +
    scale_fill_manual(values = c("#4393C3", "#D6604D"), name = "Scenario") +
    scale_colour_manual(values = c("#4393C3", "#D6604D"), guide = "none") +
    coord_cartesian(ylim = ylim) +
    labs(
      x = NULL,
      y = if (err == "abs_err")
        "Absolute error (true – estimate)"
      else
        "Relative error  (%)",
      title = "Does Choosing the Highest Log-Likelihood Help?",
      subtitle = "Grouped violins by optimizer strategy"
    ) +
    theme(
      strip.text  = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1),
      plot.background = element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      legend.background = element_blank()
    )
}

mle_accu_mmle <- plot_typ_vs_opt(mle_accu_typ, mle_accu_opt, err = "rel_err_pc", ylim = c(-200,400))
ggsave("mle_accuracy_mmle.png", mle_accu_mmle, , device = "png",
       width = 12, height = 9, dpi = "retina")

# ────────────────────────────────────────────────────────────────
# convert_differences(): long → wide  (λ, μ, K  +  errors)
# ----------------------------------------------------------------
# df          – tidy data frame returned by read_differences()
# model_type  – character tag for the “Model” column
# task_type   – character tag for the “Task”  column
# id_vars     – grouping key.  By default it is the .rds file name,
#               but you can add opt_method/integrator if you want
#               a row *per* optimiser / solver.
# ----------------------------------------------------------------
convert_differences <- function(df,
                                model_type,
                                task_type,
                                id_vars = "file") {

  stopifnot(all(c("param", "true", "rel_err_pc", "abs_err", "nnode") %in% names(df)))

  # normalise parameter names so "K" → "cap"
  df <- dplyr::mutate(df,
                      param_std = dplyr::recode(param,
                                                lambda = "lambda",
                                                mu     = "mu",
                                                K      = "cap"))

  wide <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(id_vars))) |>
    dplyr::summarise(
      lambda        = true     [param_std == "lambda"][1],
      mu            = true     [param_std == "mu"]    [1],
      cap           = true     [param_std == "cap"]   [1],
      lambda_r_diff = rel_err_pc[param_std == "lambda"][1],
      mu_r_diff     = rel_err_pc[param_std == "mu"]    [1],
      cap_r_diff    = rel_err_pc[param_std == "cap"]   [1],
      lambda_a_diff = abs_err  [param_std == "lambda"][1],
      mu_a_diff     = abs_err  [param_std == "mu"]    [1],
      cap_a_diff    = abs_err  [param_std == "cap"]   [1],
      num_nodes     = nnode[1],
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Model = model_type,
      Task  = task_type
    ) |>
    dplyr::select(lambda, mu, cap,
                  lambda_r_diff, mu_r_diff, cap_r_diff,
                  lambda_a_diff, mu_a_diff, cap_a_diff,
                  num_nodes, Model, Task)

  as.data.frame(wide)
}

library(purrr)
library(dplyr)
library(tibble)

# ─────────────────────────────────────────────────────────────────
# 1.  Loader: read one block’s differences_*.rds into a tidy tibble
# ─────────────────────────────────────────────────────────────────
load_mle_block <- function(root_path,
                           task = c("best", "typical")) {
  task   <- match.arg(task)
  subdir <- if (task == "typical")
    file.path(root_path, "DDD_MLE_TES", "NO_INIT")
  else
    file.path(root_path, "DDD_MLE_TES")

  files  <- list.files(subdir,
                       pattern = "^differences_[0-9]+\\.rds$",
                       full.names = TRUE)
  if (!length(files))
    stop("No differences_*.rds in ", subdir)

  purrr::map_dfr(files, function(f) {
    obj <- readRDS(f)

    # true parameters (always length-3 numeric)
    true_pars <- as.numeric(obj$input$truepars)

    purrr::map_dfr(obj$results, function(x) {

      # --- force everything to the right basic type ----------------
      loglik <- suppressWarnings(as.numeric(x$loglik)[1])
      est    <- suppressWarnings(as.numeric(x$est))
      rel_e  <- suppressWarnings(as.numeric(x$rel_err_pc))
      abs_e  <- suppressWarnings(as.numeric(x$abs_err))

      # guard against malformed records
      if (length(est) != 3 || anyNA(est)) {
        est   <- rep(NA_real_, 3)
        rel_e <- rep(NA_real_, 3)
        abs_e <- rep(NA_real_, 3)
      }

      tibble::tibble(
        file        = basename(f),
        start_case  = task,
        opt_method  = as.character(x$opt_method),
        integrator  = as.character(x$integrator),
        loglik      = loglik,
        param       = c("lambda", "mu", "K"),
        true        = true_pars,
        estimate    = est,
        rel_err_pc  = rel_e,
        abs_err     = abs_e,
        console_log = x$console_log,
        nnode       = as.integer(x$nnode %||% NA)
      )
    })
  })
}

# ─────────────────────────────────────────────────────────────────
# 2.  Plotter: compare estimate vs. true for the two integrators
# ─────────────────────────────────────────────────────────────────
# tag the locally best integrator for every (file × start × parameter)
tag_best_integrator <- function(df) {
  df |>
    dplyr::group_by(file, start_case, param) |>
    dplyr::mutate(best = abs(rel_err_pc) == min(abs(rel_err_pc), na.rm = TRUE)) |>
    dplyr::ungroup()
}

# ─────────────────────────────────────────────────────────────
#  Accuracy plots ─ separate scatter panels for Best & Typical
# ─────────────────────────────────────────────────────────────
plot_param_accuracy <- function(df, rel_cut = 2000, zoom_y = c(-400, 400)) {

  library(ggplot2)
  library(patchwork)

  df <- df |>
    dplyr::filter(abs(rel_err_pc) <= rel_cut) |>
    tag_best_integrator()

  # nice facet labels
  greek_lab <- ggplot2::as_labeller(
    c(lambda = "lambda",
      mu     = "mu",
      K      = "italic(K)"),
    default = ggplot2::label_parsed)

  base_scatter <- function(.data, title_suffix) {
    ggplot(.data,
           aes(true, estimate,
               colour = integrator,
               shape  = integrator,
               size   = best)) +
      geom_line(aes(group = interaction(file, param)),
                colour = "grey70", linewidth = 0.3) +
      geom_point(alpha = 0.8, stroke = 0.4) +
      geom_abline(slope = 1, intercept = 0, linewidth = 1,
                  linetype = "dashed", colour = "black") +
      scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.5),
                        guide = "none") +
      scale_colour_brewer(palette = "Dark2", name = "Integrator") +
      scale_shape_manual(values = c(16, 17), name = "Integrator") +
      facet_wrap(~ param, ncol = 3,
                 labeller = labeller(param = greek_lab),
                 scales = "free") +
      labs(x = "True value",
           y = "Estimate",
           title = paste("Estimate vs. true —", title_suffix)) +
      theme(plot.background = element_blank(), panel.background = element_blank())
  }

  p_scatter_best <- base_scatter(dplyr::filter(df, start_case == "best"),
                                 "Best start")
  p_scatter_typ  <- base_scatter(dplyr::filter(df, start_case == "typical"),
                                 "Typical start")

  p_violin <- ggplot(df,
                     aes(x = integrator,
                         y = rel_err_pc,
                         fill = start_case)) +
    geom_violin(trim = TRUE,
                alpha = 0.6,
                colour = "grey30",
                position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.12,
                 outlier.shape = NA,
                 linewidth = 0.3,
                 position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ param, ncol = 3,
               labeller = labeller(param = greek_lab),
               scales = "free") +
    coord_cartesian(ylim = zoom_y) +
    scale_fill_brewer(palette = "Set2", name = "Start case") +
    labs(x = "Integrator",
         y = "Relative error  (%)",
         title = "Relative-error distribution (≤ ±2000 %)") +
    theme(plot.background = element_blank(), panel.background = element_blank())

  # vertically stack the three panels
  p_scatter_best / p_scatter_typ / p_violin + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

# ────────────────────────────────────────────────────────────────
# 3 ▸  Example
# ────────────────────────────────────────────────────────────────
root <- "C:\\Users\\tianj\\OneDrive\\My\\Data\\COMP3"

df_best    <- load_mle_block(root, task = "best")
df_typical <- load_mle_block(root, task = "typical")

df_all <- dplyr::bind_rows(df_best, df_typical)

p_accu_int <- plot_param_accuracy(df_all)
ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\param_accurary_integrator2.png", p_accu_int, device = "png",
       width = 20, height = 28, dpi = "retina")

interactive_extract <- function(root_path) {
  p_idx <- utils::menu(c("lambda", "mu", "K"), title = "Choose parameter:")
  if (p_idx == 0) return(invisible(NULL))
  param <- c("lambda", "mu", "K")[p_idx]

  int_menu <- c("analytical", "odeint::runge_kutta_cash_karp54")
  w_idx <- utils::menu(int_menu, title = "Select winning integrator:")
  if (w_idx == 0) return(invisible(NULL))
  winner_int <- int_menu[w_idx]

  read_block <- function(subdir, start_case) {
    files <- list.files(subdir, "^differences_[0-9]+\\.rds$", full.names = TRUE)

    purrr::map_dfr(files, function(f) {
      obj   <- readRDS(f)
      res   <- purrr::keep(obj$results,
                           ~ .x$integrator %in% int_menu)
      if (length(res) != 2) return(NULL)

      idx <- as.integer(sub("^differences_(\\d+)\\.rds$", "\\1",
                            basename(f)))

      take <- function(x, field) {
        as.numeric(x[[field]][ match(param, c("lambda", "mu", "K")) ])
      }

      tibble::tibble(
        index       = idx,
        start_case  = start_case,
        integrator  = purrr::map_chr(res, "integrator"),
        estimate    = purrr::map_dbl(res, take, "est"),
        rel_err_pc  = purrr::map_dbl(res, take, "rel_err_pc"),
        abs_err     = purrr::map_dbl(res, take, "abs_err"),
        true        = obj$input$truepars[ match(param, c("lambda", "mu", "K")) ],
        brts        = list(obj$input$brts),
        initpars    = list(obj$input$initpars),
        console_log = purrr::map_chr(res, "console_log")
      )
    })
  }

  best_dir <- file.path(root_path, "DDD_MLE_TES")
  typ_dir  <- file.path(root_path, "DDD_MLE_TES", "NO_INIT")

  df_long <- dplyr::bind_rows(
    read_block(best_dir, "Best"),
    read_block(typ_dir,  "Typical")
  )

  if (!nrow(df_long)) {
    message("No comparable runs found.")
    return(invisible(NULL))
  }

  win_key <- df_long |>
    dplyr::group_by(index, start_case) |>
    dplyr::slice_min(abs(rel_err_pc), with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::filter(integrator == winner_int) |>
    dplyr::select(index, start_case)

  df_pairs <- dplyr::semi_join(df_long, win_key,
                               by = c("index", "start_case"))

  display_tbl <- dplyr::filter(df_pairs, integrator == winner_int) |>
    dplyr::select(index, start_case, integrator,
                  true, estimate, rel_err_pc, console_log) |>
    dplyr::arrange(start_case, index)

  cat("\nParameter:", param,
      "\nWinning integrator:", winner_int,
      "\nRows where it wins:\n\n")
  print(display_tbl, n = Inf)

  store_ans <- utils::menu(c("All", "Select rows", "None"),
                           title = "Store which rows (winner + loser)?")
  if (store_ans %in% c(0, 3)) {
    message("Nothing stored.")
    return(invisible(df_pairs))
  }

  save_rows <- df_pairs
  if (store_ans == 2) {
    idx_txt <- readline("Enter indices (as shown) to keep, e.g. 1,3,7: ")
    idx_in  <- as.integer(strsplit(idx_txt, ",")[[1]])
    key     <- dplyr::filter(display_tbl, index %in% idx_in) |>
      dplyr::select(index, start_case)
    save_rows <- dplyr::semi_join(df_pairs, key,
                                  by = c("index", "start_case"))
  }

  varname <- readline("Variable name for saved object: ")
  if (nzchar(varname)) {
    assign(varname, save_rows, envir = .GlobalEnv)
    message("Saved ", nrow(save_rows), " rows to object '", varname, "'.")
  } else {
    message("Invalid name – nothing stored.")
  }

  invisible(save_rows)
}

# ================================================================
#  visualize_selected()
# ----------------------------------------------------------------
#  sel_df   – data-frame returned by interactive_extract()
#             (two rows per {index × start_case}: one analytical,
#              one odeint::runge_kutta_cash_karp54)
#  title    – optional main title
# ================================================================
visualize_selected <- function(sel_df, title = NULL) {

  stopifnot(all(c("index", "start_case", "integrator",
                  "true", "estimate", "rel_err_pc") %in% names(sel_df)))

  param <- attr(sel_df, "param") %||% "chosen parameter"

  library(ggplot2); library(patchwork)

  # sort integrator factor for consistent colours / shapes
  sel_df$integrator <- factor(sel_df$integrator,
                              levels = c("analytical",
                                         "odeint::runge_kutta_cash_karp54"),
                              labels = c("Analytical", "Odeint RK54"))

  p_scatter <- ggplot(sel_df,
                      aes(x = true,
                          y = estimate,
                          colour = integrator,
                          shape  = integrator)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "black") +
    geom_line(aes(group = interaction(index, start_case)),
              colour = "grey70", linewidth = 0.3) +
    geom_point(size = 2.5, alpha = 0.8) +
    facet_wrap(~ start_case, ncol = 1, scales = "free") +
    scale_colour_brewer(palette = "Dark2", name = "Integrator") +
    scale_shape_manual(values = c(16, 17), name = "Integrator") +
    coord_cartesian(ylim = c(-1, 30)) +
    labs(x = "True value",
         y = "Estimate",
         title = "Estimate vs. true",
         subtitle = paste("Parameter:", param)) +
    theme_bw()

  p_violin <- ggplot(sel_df,
                     aes(x = integrator,
                         y = rel_err_pc,
                         fill = integrator)) +
    geom_violin(alpha = 0.6, colour = "grey30", width = 0.9) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ start_case, ncol = 1, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    coord_cartesian(ylim = c(-2000,2000)) +
    labs(x = NULL,
         y = "Relative error  (%)",
         title = "Relative error distribution") +
    theme_bw()

  full_title <- title %||% paste("Selected data –", param)

  (p_scatter + p_violin) +
    patchwork::plot_annotation(title = full_title) &
    theme(plot.title = element_text(hjust = 0.5))
}

root <- "path to unzipped folder" # should be "COMP", one level higher than "DDD_MLE_TES"
interactive_extract(root) # Interactive selector

#### Assume you stored selection as "my_runs" ###
#### and you choose "lambda", and "analytical" as winner ###
attr(my_runs, "param") <- "lambda"

# Generate plot
p <- visualize_selected(my_runs,
                   title = "Winner–loser comparison (λ)")

ggsave("C:\\Users\\tianj\\OneDrive\\My\\Projects\\PhD\\Thesis-Chapter2\\revision\\figure\\param_accurary_lambda_ana2.png", p, device = "png",
       width = 20, height = 28, dpi = "retina")

# True pars, brts and initial pars are also in my_runs
names(my_runs)

### Script contributed by Koen (B.)
# Extract runs with errors
bla <- lapply(which(!unlist(success)), FUN = function(z){
  y <- my_runs$console_log[z]
  tmp <- strsplit(y,split="\n")[[1]]
  data.frame(errmessage = paste(tmp[max(1,length(tmp)-2):length(tmp)],collapse="\n"),
             integrator = my_runs$integrator[z],
             case = my_runs$start_case[z],
             index=z)
})
bla2 <- do.call(rbind,bla)

# Overview of errors
table(bla2$case,bla2$errmessage)
table(bla2$integrator,bla2$errmessage)

# Here, find index of runs with typical, odeint and -Inf error message
bla2[bla2$case=='Typical' & bla2$integrator=='odeint::runge_kutta_cash_karp54' & grepl("-Inf",bla2$errmessage),]

# Show information of data of interest
# with names(my_runs)
cat(lapply(my_runs, FUN = function(x) x[1022])$integrator)

# Note that here the index refer to the file name of the RDSs
# e.g. differences_XXXX.rds
# If the below index is 100, and if in the previous step you
# extract with case == 'Best', then you can try to just load
# data = readRDS("DDD_MLE_TES/differences_100.rds)
# If with case == 'Typical', then read:
# data = readRDS("DDD_MLE_TES/NO_INIT/differences_100.rds)
cat(lapply(my_runs, FUN = function(x) x[836])$index)

# Run it
# DDD::dd_ML(brts=...,
#       initparsopt = ...,
#            idparsopt   = c(1, 2, 3),
#            btorph      = 0,
#            soc         = 2,
#            cond        = 1,
#            ddmodel     = 1,
#            num_cycles  = Inf,
#            optimmethod = "simplex",
#            methode="odeint::runge_kutta_cash_karp54")