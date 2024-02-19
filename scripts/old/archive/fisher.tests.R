rm(list=ls())


gene.type.raw <- read.csv("~/GitHub/Bifidobacterium/4_gene.types.chisq.csv")

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
gene.type$mut_type <- factor(gene.type$mut_type, levels = c("Transport", "Metabolism","Other"))
gene.type <- gene.type %>% mutate(prop2 = sprintf("%0.2f", prop))

# Reorder for plotting
order.1 <- c(2, 1, 3, 5, 4, 6, 8, 7, 9)
gene.type <- data.frame(gene.type, order.1)
gene.type <- gene.type[order(gene.type$order.1),] 

# Make table for contingency analyses
gene.type.mat <- matrix(gene.type$count, ncol = 3, byrow = T)
rownames(gene.type.mat) <- c("F", "G-F5","G-F24")
colnames(gene.type.mat) <- c("Transport", "Metabolism", "Other")
gene.type.tab <- as.table(gene.type.mat)
gene.type.tab.margins <- addmargins(gene.type.tab)

# X-squared = 3.6573, df = 4, p-value = 0.4544
gene.type.chi <- chisq.test(gene.type.tab)

#Fisher test. p = 0.5154
fisher.test(gene.type.tab)

# Posthoc analysis
gene.type.z <- as.data.frame(gene.type.chi$stdres)
gene.type.x2 <- gene.type.z$Freq^2
gene.type.p <- pchisq(gene.type.x2, df = 1, lower.tail = FALSE)
gene.type.p.adj <- p.adjust(gene.type.p, method="BH")
gene.type.post.hoc <- data.frame(gene.type.z, gene.type.x2, gene.type.p, gene.type.p.adj)
colnames(gene.type.post.hoc) <- c("type", "strain", "z", "chi2", "p", "p.adj")

#type     strain           z        chi2          p     p.adj
#1     F  Transport  0.74717513 0.558270677 0.45495785 0.6824368
#2  G-F5  Transport  0.08574213 0.007351713 0.93167141 0.9316714
#3 G-F24  Transport -0.85742129 0.735171261 0.39121210 0.6824368
#4     F Metabolism -1.40685881 1.979251701 0.15946926 0.6824368
#5  G-F5 Metabolism -0.31586903 0.099773243 0.75210192 0.9147846
#6 G-F24 Metabolism  1.76886655 3.128888889 0.07691615 0.6824368
#7     F      Other  0.80104099 0.641666667 0.42310792 0.6824368
#8  G-F5      Other  0.23637474 0.055873016 0.81314190 0.9147846
#9 G-F24      Other -1.06368631 1.131428571 0.28747083 0.6824368



##################################################################

##################################################################
#Ignoring irrelevant genes, only looking at metabolism vs. transport