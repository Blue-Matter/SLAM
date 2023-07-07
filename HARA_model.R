library(SLAM)
library(ggplot2)

Simulation <- Simulate(LifeHistory_Pulse, Exploitation, nsim=2)

sim <- 1
n_months_total <- ncol(Simulation$Exploitation$Sel_at_Age[sim,,])
sel_at_age <- Simulation$Exploitation$Sel_at_Age[sim,,n_months_total]

utilpow_vec <- seq(0.2, 1, by=0.2)
util_list <- list()

for (i in seq_along(utilpow_vec)) {
  OM <- calculate_optimal_fishing(R0_m=Simulation$LifeHistory$R0_m,
                                  steepness=Simulation$LifeHistory$steepness,
                                  Weight_Age_Mean=Simulation$LifeHistory$Weight_Age_Mean,
                                  Maturity_at_Age=Simulation$LifeHistory$Maturity_at_Age,
                                  M_at_Age=Simulation$LifeHistory$M_at_Age,
                                  Post_Spawning_Mortality=Simulation$LifeHistory$Post_Spawning_Mortality,
                                  sel_at_age=sel_at_age,
                                  opt_type=1,
                                  utilpow=utilpow_vec[i])

  util_list[[i]] <- data.frame(utilpow=utilpow_vec[i], Month=month.abb, Catch=OM$predCB, Effort=OM$F_m)
}



df <- data.frame(Month=month.abb, Rec=Simulation$LifeHistory$R0_m)
df$Month <- factor(df$Month, levels=month.abb, ordered = TRUE)


ggplot(df, aes(x=Month, y=Rec, group=1)) +
  geom_line()  +
  expand_limits(y=0) +
  labs(y='Relative Recruitment') +
  theme_bw()

ggsave('img/HARA/Recruitment.png', width=4, height=3)


DF <- do.call('rbind', util_list)
DF$Month <- factor(DF$Month, levels=month.abb, ordered = TRUE)


txt_df <- DF %>% group_by(utilpow) %>% summarize(Catch=sum(Catch))
txt_df$Catch <- txt_df$Catch/max(txt_df$Catch)
txt_df$Catch <- round(txt_df$Catch,2)
txt_df$name <- 'Catch'
#
# ggplot(DF, aes(x=Month, y=Catch, group=utilpow)) +
#   facet_wrap(~utilpow) +
#   expand_limits(y=0) +
#   geom_line() +
#   geom_text(data=txt_df, aes(x='Jan', y=Inf, label=Catch),hjust = 0, vjust = 1) +
#   theme_bw() +
#   labs(y='Relative Catch')
#
# ggplot(DF, aes(x=Month, y=Effort, group=utilpow)) +
#   facet_wrap(~utilpow) +
#   expand_limits(y=0) +
#   geom_line() +
#   geom_text(data=txt_df, aes(x='Jan', y=Inf, label=Catch),hjust = 0, vjust = 1) +
#   theme_bw() +
#   labs(y='Relative Effort')

DF2 <- DF %>% tidyr::pivot_longer(., cols=3:4)

ggplot(DF2, aes(x=Month, y=value, group=utilpow)) +
  facet_grid(name~utilpow, scales='free') +
  expand_limits(y=0) +
  geom_line() +
  geom_text(data=txt_df, aes(x='Jan', y=Inf, label=Catch),hjust = 0, vjust = 1) +
  theme_bw() +
  labs(y='Relative Value') +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

ggsave('img/HARA/HARA.png', width=8, height=4)
