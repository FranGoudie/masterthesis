#### File to look at number of iterations needed in bootstrap and in Sim Study ##
## Also has code for all graphs used in master thesis ##

## Load libraries ##
library(tidyverse)
library(ggtext)
library(patchwork)

# Do this the same but in parrallel so that it is faster
myCluster <- parallel::makeCluster(detectCores(), type = "PSOCK",
                                   output = "")
registerDoParallel(myCluster)
clusterExport(myCluster, varlist = c("X", 'Y', 'Z',
                                     'hot_deck', 'quality_SM',
                                     "A", "B", "C", "CinB", "CinA"))
k <- 1
l <- 4
iterations <- list()
num <- c(10,20,30,40,50,60,70,80,90,100,130,160,190,220)
myCluster <- parallel::makeCluster(detectCores(), type = "PSOCK",
                                   output = "")
registerDoParallel(myCluster)
clusterExport(myCluster, varlist = c("X", 'Y', 'Z',
                                     'dist_hotdeck', 'quality_SM',
                                     "A", "B", "C", "CinB", "CinA"))

iterations_to <- foreach(i=num,
                         .packages = c("tidyverse", "StatMatch"))%dopar%{
                           quality_est <- quality_SM(A = A[-CinA, c(Y,X)], 
                                                     B = B[-CinB, c(Z,X)],
                                                     C = C[,c(Z,Y,X)],
                                                     SM = hot_deck,
                                                     Y = Y, Z = Z, X = X,
                                                     R = i, aux = 1)
                           n_iter_bias_2 <- quality_est[[1]]
                           n_iter_var_2 <- quality_est[[3]]
                           
                           iter <- list(n_iter_bias_2, n_iter_var_2)
                           
                           return(iter)}
stopCluster(myCluster)

for(i in 1:12){
  n_iter_bias_2[i] <- iterations_to[[i]][1]
  n_iter_var_2[i] <- iterations_to[[i]][2]
}

# Look at each cell of contingency table
itg_bias <- list()
itg_var <- list()
for(i in 1:15){
  itg_bias[[i]] <- rep(NA, times = 12)
  itg_var[[i]] <- rep(NA, times = 12)
  for(j in 1:12){
    itg_var[[i]][j] <- unlist(n_iter_var_2)[(15*(j-1)+i)]
    itg_bias[[i]][j] <- unlist(n_iter_bias_2)[(15*(j-1)+i)]
  }
}

# iterations
iter <- c(10,20,30,40,50,60,70,80,90,100,130,160,190,220)

plot_bv <- function(x){plot(iter, x)}

varplots_2 <- lapply(itg_var, plot_bv)
biasplots_2 <- lapply(itg_bias, plot_bv)

# 4 plot showing the iterations of number 1,4,7,15
plot_data <- data.frame(iter = iter,
                        bias_15 = itg_bias[[15]],
                        bias_7 = itg_bias[[7]],
                        bias_10 = itg_bias[[10]],
                        bias_1 = itg_bias[[1]])

p1 <- ggplot(plot_data, aes(x = iter, y = bias_15)) +
  geom_point() + geom_line() + ylim(0.001, 0.003) +
  labs(x = "Number of iterations",
       y = "Bias of Cell (5,3)")
p2 <- ggplot(plot_data, aes(x = iter, y = bias_10)) +
  geom_point() + geom_line() + ylim(-0.012, -0.007) +
  labs(x = "Number of iterations",
       y = "Bias of Cell (5,2)")
p3 <- ggplot(plot_data, aes(x = iter, y = bias_7)) +
  geom_point() + geom_line() + ylim(0.016, 0.023) +
  labs(x = "Number of iterations",
       y = "Bias of Cell (2,2)")
p4 <- ggplot(plot_data, aes(x = iter, y = bias_1)) +
  geom_point() + geom_line() + ylim(0.002,0.007) +
  labs(x = "Number of iterations",
       y = "Bias of Cell (1,1)")
p1+p2+p3+p4 + plot_annotation(
  title = "Bias over Bootstrap iterations"
)

### Simulation iterations ###
# Use quality_data_dist
load(paste0(sourcepath, "quality_data_dist.Rdata"))
# Uses t from the simulation_study.R file based on the true contigency table
# of Gez_sim data

#Look at iterations of true_bias
iter_true_bias <- rep(list(rep(NA, times = 100)), 15)
for(i in 1:100){
  true_bias_list_it <- lapply(est_dist[1:i], FUN = function(x){x - t})
  true_bias_red_it <- Reduce('+', true_bias_list_it)
  true_bias_it <- (1/i)*true_bias_red_it
  for(j in 1:15){
    iter_true_bias[[j]][i] <- as.vector(true_bias_it)[j]
  }
}

plot_iter <- function(x){plot(1:100, x, pch = 20)}

true_plots <- lapply(iter_true_bias, plot_iter)

iter_true_var <- rep(list(rep(NA, times = 100)), 15)
for(i in 1:100){
  mean_est <- (1/i)*Reduce('+', est_dist[1:i])
  true_var_list_it <- lapply(est_dist[1:i], FUN = function(x){(x - mean_est)^2})
  true_var_red_it <- Reduce('+', true_var_list_it)
  true_var_it <- (1/(i-1))*true_var_red_it
  for(j in 1:15){
    iter_true_var[[j]][i] <- as.vector(true_var_it)[j]
  }
}

plot_iter <- function(x){plot(1:100, x, pch = 20)}

true_plots <- lapply(iter_true_var, plot_iter)

it_data <- data.frame(iter = 1:100,
                      iter_bias_8 = iter_true_bias[[8]],
                      iter_bias_15 = iter_true_bias[[15]],
                      iter_var_6 = iter_true_var[[6]],
                      iter_var_11 = iter_true_var[[11]],
                      iter_var_12 = iter_true_var[[12]],
                      iter_var_14 = iter_true_var[[14]])

i1 <- ggplot(it_data, aes(x = iter, y = iter_bias_8)) +
  geom_point() + ylim(-0.001, 0.004) +
  labs(x = "Number of iterations",
       y = "True Bias of Cell (3,3)")
i2 <- ggplot(it_data, aes(x = iter, y = iter_bias_15)) +
  geom_point() + ylim(0, 0.003) +
  labs(x = "Number of iterations",
       y = "True Bias of Cell (5,3)")
i3 <- ggplot(it_data, aes(x = iter, y = iter_var_14)) +
  geom_point() + ylim(6.2e-07, 2.3e-06) +
  labs(x = "Number of iterations",
       y = "True Variance of Cell (4,3)")
i4 <- ggplot(it_data, aes(x = iter, y = iter_var_11)) +
  geom_point() + ylim(0, 8e-07) +
  labs(x = "Number of iterations",
       y = "True Variance of Cell (1,3)")
i1 + i2 + i3 + i4 + plot_annotation(
  title = "True bias and variance over Simulation iterations"
)

###################################
### Plots used in master thesis ###
###################################

# Uses results from simulation study

# Histograms of bias for random hotdeck procedure
hist_bias_ran <- rep(list(rep(NA, times = 100)), 15)
for(i in 1:100){
  for(j in 1:15){
    hist_bias_ran[[j]][i] <- unlist(quality_est_bias_random)[(15*(i-1)+j)]
  }
}

hist_bias_ran <- as.data.frame(hist_bias_ran)
colnames(hist_bias_ran)<- letters[1:15]
hist_bias_ran <- gather(hist_bias_ran)

bias_ran_mean <- hist_bias_ran %>% group_by(key) %>%
  summarise(mean = mean(value))

hist_bias_ran <- full_join(hist_bias_ran, bias_ran_mean, by = "key")

bias_true_ran_plot <- data.frame(key = letters[1:15],
                                 true_bias = as.vector(true_bias_ran)) 

hist_bias_ran <- full_join(hist_bias_ran, bias_true_ran_plot, by = "key")

names <- c(
  'a' = 'Cell (1,1)',
  'b' = 'Cell (2,1)', 
  'c' = 'Cell (3,1)',
  'd' = 'Cell (4,1)',
  'e' = 'Cell (5,1)',
  'f' = 'Cell (1,2)',
  'g' = 'Cell (2,2)',
  'h' = 'Cell (3,2)',
  'i' = 'Cell (4,2)',
  'j' = 'Cell (5,2)',
  'k' = 'Cell (1,3)',
  'l' = 'Cell (2,3)',
  'm' = 'Cell (3,3)',
  'n' = 'Cell (4,3)',
  'o' = 'Cell (5,3)'
)

ggplot(data = hist_bias_ran) +
  geom_histogram(aes(x = value), bins = 20, color = "black", fill = 'grey') +
  facet_wrap(~key, scales = "free_x", ncol = 3, 
             labeller = as_labeller(names), dir = 'v') +
  geom_vline(aes(xintercept=mean, colour = "black", linetype = "dashed")) +
  geom_vline(aes(xintercept=true_bias, colour = "red", linetype = "solid")) + 
  labs(y = "Count",
       x = "Bias estimate",
       title = "Histogram bootstrap
       bias estimates over the simulation
       using the random hotdeck procedure",
       linetype = "Lines",
       colour = "Lines") +
  scale_colour_manual(name = "Lines",
                      breaks = c("black", "red"),
                      labels = c("Mean of Bias estimates",
                                 "True Bias"),
                      values = c(black = "black",
                                 red = "red")) +
  scale_linetype_manual(name = "Lines",
                        values = c(1, 2),
                        labels = c("Mean of Bias estimates",
                                   "True Bias")) +
  theme(plot.title = element_textbox_simple(height = 0.1),
        panel.spacing.x = unit(6, "mm"))

# Histograms of relative bias for random hotdeck procedure
hist_rel_bias_ran <- rep(list(rep(NA, times = 100)), 15)
for(i in 1:100){
  for(j in 1:15){
    hist_rel_bias_ran[[j]][i] <- unlist(quality_est_rel_bias_random)[(15*(i-1)+j)]
  }
}

hist_rel_bias_ran <- as.data.frame(hist_rel_bias_ran)
colnames(hist_rel_bias_ran)<- letters[1:15]
hist_rel_bias_ran <- gather(hist_rel_bias_ran)

hist_rel_bias_ran <- hist_rel_bias_ran[is.finite(hist_rel_bias_ran$value),]

bias_rel_ran_mean <- hist_rel_bias_ran %>%
  group_by(key) %>%
  summarise(mean = mean(value))

hist_rel_bias_ran <- full_join(hist_rel_bias_ran, bias_rel_ran_mean, by = "key")

bias_rel_true_ran_plot <- data.frame(key = letters[1:15],
                                     true_bias = as.vector(true_rel_bias_ran))

hist_rel_bias_ran <- full_join(hist_rel_bias_ran, bias_rel_true_ran_plot,
                               by = "key")

ggplot(data = hist_rel_bias_ran) +
  geom_histogram(aes(x = value), bins = 20, color = "black", fill = 'grey') +
  facet_wrap(~key, scales = "free_x", ncol = 3, 
             labeller = as_labeller(names), dir = 'v') +
  geom_vline(aes(xintercept=mean, colour = "black", linetype = "dashed")) +
  geom_vline(aes(xintercept=true_bias, colour = "red", linetype = "solid")) + 
  labs(y = "Count",
       x = "Bias estimate",
       title = "Histogram bootstrap
       relative bias estimates over the simulation
       using the random hotdeck procedure",
       linetype = "Lines",
       colour = "Lines") +
  scale_colour_manual(name = "Lines",
                      breaks = c("black", "red"),
                      labels = c("Mean of relative bias estimates",
                                 "True Relative Bias"),
                      values = c(black = "black",
                                 red = "red")) +
  scale_linetype_manual(name = "Lines",
                        values = c(1, 2),
                        labels = c("Mean of relative bias estimates",
                                   "True Relative Bias")) +
  theme(plot.title = element_textbox_simple(height = 0.15),
        panel.spacing.x = unit(6, "mm"))

# Histograms of true estimation error for distance hot deck procedure 
hist_tru_est_dist <- rep(list(rep(NA, times = 100)), 15)
for(i in 1:100){
  for(j in 1:15){
    hist_tru_est_dist[[j]][i] <- unlist(true_est_err_list_dist)[(15*(i-1)+j)]
  }
}

hist_tru_est_dist <- as.data.frame(hist_tru_est_dist)
colnames(hist_tru_est_dist)<- letters[1:15]
hist_tru_est_dist <- gather(hist_tru_est_dist)

est_dist_mean <- hist_tru_est_dist %>%
  group_by(key) %>%
  summarise(mean = mean(value))

hist_tru_est_dist <- full_join(hist_tru_est_dist, est_dist_mean, by = "key")

true_bias_plot_dist <- data.frame(key = letters[1:15],
                                  true_bias = as.vector(true_bias_dist))

hist_tru_est_dist <- full_join(hist_tru_est_dist, true_bias_plot_dist,
                               by = "key")

ggplot(data = hist_tru_est_dist) +
  geom_histogram(aes(x = value), bins = 20, color = "black", fill = 'grey') +
  facet_wrap(~key, scales = "free_x", ncol = 3, 
             labeller = as_labeller(names), dir = 'v') +
  geom_vline(aes(xintercept=mean, colour = "black", linetype = "dashed")) +
  geom_vline(aes(xintercept=true_bias, colour = "red", linetype = "solid")) + 
  labs(y = "Count",
       x = "Estimation errror",
       title = "Cramer's V for income and proxy income with distance hot deck",
       linetype = "Lines",
       colour = "Lines") +
  scale_colour_manual(name = "Lines",
                      breaks = c("black", "red"),
                      labels = c("Average true estimation error",
                                 "True Bias"),
                      values = c(black = "black",
                                 red = "red")) +
  scale_linetype_manual(name = "Lines",
                        values = c(1, 2),
                        labels = c("Average true estimation error",
                                   "True Bias")) +
  theme(plot.title = element_textbox_simple(height = 0.15),
        panel.spacing.x = unit(6, "mm"))


# Histograms of true relative estimation error for distance hot deck procedure 
hist_tru_rel_dist <- rep(list(rep(NA, times = 100)), 15)
for(i in 1:100){
  for(j in 1:15){
    hist_tru_rel_dist[[j]][i] <- unlist(true_rel_err_list_dist)[(15*(i-1)+j)]
  }
}

hist_tru_rel_dist <- as.data.frame(hist_tru_rel_dist)
colnames(hist_tru_rel_dist)<- letters[1:15]
hist_tru_rel_dist <- gather(hist_tru_rel_dist)

rel_dist_mean <- hist_tru_rel_dist %>%
  group_by(key) %>%
  summarise(mean = mean(value))

hist_tru_rel_dist <- full_join(hist_tru_rel_dist, rel_dist_mean, by = "key")

true_rel_bias_plot_dist <- data.frame(key = letters[1:15],
                                      true_bias = as.vector(true_rel_bias_dist))

hist_tru_rel_dist <- full_join(hist_tru_rel_dist, true_rel_bias_plot_dist,
                               by = "key")

ggplot(data = hist_tru_rel_dist) +
  geom_histogram(aes(x = value), bins = 20, color = "black", fill = 'grey') +
  facet_wrap(~key, scales = "free_x", ncol = 3, 
             labeller = as_labeller(names), dir = 'v') +
  geom_vline(aes(xintercept=mean, colour = "black", linetype = "dashed")) +
  geom_vline(aes(xintercept=true_bias, colour = "red", linetype = "solid")) + 
  labs(y = "Count",
       x = "Relative estimation error",
       title = "Cramer's V for income and proxy income with distance hot deck",
       linetype = "Lines",
       colour = "Lines") +
  scale_colour_manual(name = "Lines",
                      breaks = c("black", "red"),
                      labels = c("Average true relative estimation error",
                                 "True Relative Bias"),
                      values = c(black = "black",
                                 red = "red")) +
  scale_linetype_manual(name = "Lines",
                        values = c(1, 2),
                        labels = c("Average true relative estimation error",
                                   "True Relative Bias")) +
  theme(plot.title = element_textbox_simple(height = 0.15),
        panel.spacing.x = unit(6, "mm"))


