library(GOSemSim)
library(GO.db)
library(org.Hs.eg.db)
library(org.At.tair.db)

BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.At.tair.db")

set.seed(123)
ontology <- "BP"
n_terms <- 10

# Obtener términos GO
all_go <- keys(GO.db, keytype = "GOID")
go_terms <- all_go[sapply(all_go, function(go) {
  ont <- Ontology(go)
  !is.na(ont) && ont == ontology
})]

go_sample <- sample(go_terms, n_terms)

sem_go <- godata(
  ont = ontology,
  computeIC = FALSE
)

sim_go <- mgoSim(
  go_sample,
  go_sample,
  semData = sem_go,
  measure = "Wang",
  combine = NULL
)

# convertir a data.frame
df_go <- as.data.frame(sim_go)

# guardar
write.csv(df_go, "wang_similarity_GO_terms.csv")

cat("GO terms Wang matrix guardada\n")


# seleccionar 10 genes humanos al azar (Entrez ID)
human_symbols <- c(
  "VSTM2L", "N4BP1", "GIT2", "OMA1", "AP2A1",
  "SSU72", "KLHL18", "FOXH1", "MRPS15", "TTC19"
)

human_entrez <- mapIds(
  org.Hs.eg.db,
  keys = human_symbols,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# eliminar NA (muy importante)
human_entrez <- human_entrez[!is.na(human_entrez)]

sem_human <- godata(
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  computeIC = FALSE
)

sim_human <- mgeneSim(
  genes = human_entrez,
  semData = sem_human,
  measure = "Wang",
  combine = "BMA"
)

df_human <- as.data.frame(sim_human)

rownames(df_human) <- names(human_entrez)
colnames(df_human) <- names(human_entrez)

write.csv(df_human, "wang_similarity_human_genes.csv")
cat("✔ Human genes Wang matrix guardada con SYMBOL\n")


# seleccionar 10 genes AGI al azar
arab_genes <- c(
  'AT4G16610', 'ELI3-1', 'CAM2', 'CCR1', 'TSA1',
  'ECS1', 'HOG1', 'AT-P4H-1', 'HLECRK', 'LTI30'
)

arab_entrez <- mapIds(
  org.At.tair.db,
  keys = arab_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

arab_entrez <- arab_entrez[!is.na(arab_entrez)]

sem_arab <- godata(
  OrgDb = org.At.tair.db,
  ont = "BP",
  computeIC = FALSE
)

sim_arab <- mgeneSim(
  genes = arab_entrez,
  semData = sem_arab,
  measure = "Wang",
  combine = "BMA"
)

df_arab <- as.data.frame(sim_arab)

rownames(df_arab) <- names(arab_entrez)
colnames(df_arab) <- names(arab_entrez)

write.csv(df_arab, "wang_similarity_arabidopsis_genes.csv")
cat("✔ Arabidopsis genes Wang matrix guardada con AGI\n")