require("png")
require("dplyr")
require("vioplot")
require("grid")
require("dplyr")
require("knitr")
require("extrafont")

rm(list=ls())


gene.type.raw <- read.csv("~/GitHub/Bifidobacterium/4_gene.types.chisq_2.categs.of.interest.csv")

# Preserve row names
gene.type <- as_tibble(gene.type.raw)

# Set data as factors and rename levels
gene.type$Treatment <- as.factor(gene.type$Treatment)
gene.type$mut_type <- as.factor(gene.type$mut_type)
levels(gene.type$Treatment)[levels(gene.type$Treatment)=="F"] <- "F"
levels(gene.type$Treatment)[levels(gene.type$Treatment)=="P"] <- "G-F5"
levels(gene.type$Treatment)[levels(gene.type$Treatment)=="X"] <- "G-F24"

levels(gene.type$mut_type)[levels(gene.type$mut_type)=="Carbohydrate metabolism"] <- "Metabolism"
levels(gene.type$mut_type)[levels(gene.type$mut_type)=="Carbohydrate transport"] <- "Transport"

gene.type$Treatment <- factor(gene.type$Treatment, levels=c("F","G-F5","G-F24"))
gene.type$mut_type <- factor(gene.type$mut_type, levels = c("Transport", "Metabolism"))
gene.type <- gene.type %>% mutate(prop2 = sprintf("%0.2f", prop))

# Reorder for plotting
#order.1 <- c(2, 1, 3, 5, 4, 6, 8, 7, 9)
#gene.type <- data.frame(gene.type, order.1)
gene.type <- gene.type[order(gene.type$order),] 

# Make table for contingency analyses
gene.type.mat <- matrix(gene.type$count, ncol = 2, byrow = T)
rownames(gene.type.mat) <- c("F", "G-F5","G-F24")
colnames(gene.type.mat) <- c("Transport", "Metabolism")
gene.type.tab <- as.table(gene.type.mat)
gene.type.tab.margins <- addmargins(gene.type.tab)

# X-squared = 2.3413, df = 2, p-value = 0.3102
gene.type.chi <- chisq.test(gene.type.tab)

#Fisher test. p = 0.2421
fisher.test(gene.type.tab)

# Posthoc analysis
gene.type.z <- as.data.frame(gene.type.chi$stdres)
gene.type.x2 <- gene.type.z$Freq^2
gene.type.p <- pchisq(gene.type.x2, df = 1, lower.tail = FALSE)
gene.type.p.adj <- p.adjust(gene.type.p, method="BH")
gene.type.post.hoc <- data.frame(gene.type.z, gene.type.x2, gene.type.p, gene.type.p.adj)
colnames(gene.type.post.hoc) <- c("type", "strain", "z", "chi2", "p", "p.adj")

#   type     strain          z       chi2         p     p.adj
#1     F  Transport  1.2777531 1.63265306 0.2013365 0.3020047
#2  G-F5  Transport  0.2129589 0.04535147 0.8313591 0.8313591
#3 G-F24  Transport -1.3944334 1.94444444 0.1631868 0.3020047
#4     F Metabolism -1.2777531 1.63265306 0.2013365 0.3020047
#5  G-F5 Metabolism -0.2129589 0.04535147 0.8313591 0.8313591
#6 G-F24 Metabolism  1.3944334 1.94444444 0.1631868 0.3020047
