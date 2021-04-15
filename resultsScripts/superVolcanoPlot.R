library(tidyverse)
minOdds = 2
minP = .0000137

superZ = read_csv("superPropTestVolc.csv")
superOdds = read_csv("superEffectForVolc.csv")
superData = inner_join(superZ, superOdds, by=c("SNP", "risk allele", "population1", "population2"))

ggplot(data=superData, aes(x=superData$`odds ratio`, y=superData$p_value)) + geom_point()

newData = mutate(superData, logOdds = log2(`odds ratio`))
de = newData
# add a column of NAs
de$diffexpressed <- "Low Odds"
de$diffexpressed[de$logOdds > minOdds & de$p_value < minP] <- "High Odds"
de$diffexpressed[de$logOdds < -minOdds & de$p_value < minP] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=logOdds, y=-log10(p_value), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=minOdds, col="red") +
  geom_hline(yintercept=-log10(minP), col="red")
print(p2)

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "High Odds", "Low Odds")
p3 <- p2 + scale_colour_manual(values = mycolors)

de$delabel <- NA
de$delabel[de$diffexpressed != "Low Odds"] <- de$SNP[de$diffexpressed != "Low Odds"]

ggplot(data=de, aes(x=logOdds, y=-log10(p_value), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=de, aes(x=logOdds, y=-log10(p_value), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  xlab("log2(odds ratio)") + 
  ylab("-log10(p-value)") +
  ggtitle("Superpopulation Volcano Plot") +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=minOdds, col="red") +
  geom_hline(yintercept=-log10(minP), col="red")

