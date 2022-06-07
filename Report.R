Len_Age <- data.frame(Ages=Ages,
                      Mean=LatAge,
                      Upper=LenSD*1.96+LatAge,
                      Lower=LatAge-LenSD*1.96)

library(ggplot2)
ggplot(Len_Age, aes(x=Ages, y=Mean)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), color='lightgray') +
  geom_line(size=2) +
  theme_classic() +
  labs(x='Age (month)',
       y='Mantle Length (mm)')

List <- list()
for (a in 1:nages) {
  age <- a - 1
  Age <- age
  Length <- Len_Mids
  Prob <- ALK[a,]
  List[[a]] <- data.frame(Age=Age, Length=Length, Prob=Prob)
}
df <- do.call('rbind', List)
df$Age <- factor(df$Age)
ggplot(df, aes(x=Length, y=Prob, color=Age)) +
  geom_line()   +
  theme_classic() +
  labs(x='Mantle Length (mm)',
       y='Probability',
       color='Age (month)')
