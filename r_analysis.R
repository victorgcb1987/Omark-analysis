library(lme4)

df <- read.table("hog_species_for_GLMM.tsv", header=TRUE, sep="\t")

df$logTE <- log1p(df$TE_score)
df$logGene <- log1p(df$Gene_score)
df$logPseudo <- log1p(df$Pseudo_score)

# modelo mixto logístico
m <- glmer(AnyInconsistent ~ logTE + logGene + logPseudo + (1|HOG_ID),
           data=df, family=binomial)

summary(m)