# Librerías - En esta sección se asegura que el script posea los paquetes
# necesarios para funcionar.
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

ensure_pkg("GOSemSim", bioc = TRUE)
ensure_pkg("GO.db", bioc = TRUE)
ensure_pkg("org.Hs.eg.db", bioc = TRUE)
ensure_pkg("org.At.tair.db", bioc = TRUE)

# Parámetros - Aquí se setea una semilla para replicabilidad de selección de términos,
# la ontología que se desea probar y la cantidad de términos de prueba. Asimismo, se 
# establece el conjunto de genes de pruebas para las siguientes secciones.
# Los genes se seleccionaron de forma arbitraria teniendo en cuenta 
set.seed(123)
ontology = "BP"
output_dir <- "C:/Users/benja/Desktop/workspace/GO3_Study/R_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
n_terms <- 10

human_symbols <- c(
  "VSTM2L", "N4BP1", "GIT2", "OMA1", "AP2A1",
  "SSU72", "KLHL18", "FOXH1", "MRPS15", "TTC19"
)

arab_genes <- c(
  'AT4G16610', 'ELI3-1', 'CAM2', 'CCR1', 'TSA1',
  'ECS1', 'HOG1', 'AT-P4H-1', 'HLECRK', 'LTI30'
)

# GOTerms - Se obtienen los GOTerms de la estructura asociada 
# a la estructura disponible por la versión (R 4.4.2).
all_go <- keys(GO.db, keytype = "GOID")
go_terms <- all_go[sapply(all_go, function(go) {
  ont <- Ontology(go)
  !is.na(ont) && ont == ontology
})]
go_sample <- sample(go_terms, n_terms)

####################################### PRUEBA 1 - Cálculo de similitud semántica entre GOTerms elegidos.
# Datos de los términos.
sem_go <- godata(
  ont = ontology,
  computeIC = FALSE
)

# Matrices de similitud.
sim_go_wang <- mgoSim(
  go_sample,
  go_sample,
  semData = sem_go,
  measure = "Wang",
  combine = NULL
)


# convertir a data.frame
df_go_wang <- as.data.frame(sim_go_wang)

# guardar
write.csv(df_go_wang, file.path(output_dir, "wang_similarity_GO_terms.csv"))

# Mensaje.
cat("GO terms Wang matrix guardada\n")
####################################### PRUEBA 1 - Cálculo de similitud semántica entre GOTerms elegidos.

####################################### PRUEBA 2 - Cálculo de similitud semántica entre genes del homo sapiens.
human_entrez <- mapIds(
  org.Hs.eg.db,
  keys = human_symbols,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Eliminar NA (muy importante)
human_entrez <- human_entrez[!is.na(human_entrez)]

# Datos de los términos.
sem_human <- godata(
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  computeIC = FALSE
)

# Matrices de similitud.
sim_human_wang <- mgeneSim(
  genes = human_entrez,
  semData = sem_human,
  measure = "Wang",
  combine = "BMA"
)

# convertir a data.frame
df_human_wang <- as.data.frame(sim_human_wang)

# Se mantienen los nombres de los símbolos.
rownames(df_human_wang) <- names(human_entrez)
colnames(df_human_wang) <- names(human_entrez)

# guardar y mensaje.
write.csv(df_human_wang, file.path(output_dir, "wang_similarity_human_genes.csv"))
cat("✔ Human genes Wang matrix guardada con SYMBOL\n")
####################################### PRUEBA 2 - Cálculo de similitud semántica entre genes del homo sapiens.

####################################### PRUEBA 3 - Cálculo de similitud semántica entre genes de la arabdosis Thaliana.

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

sim_arab_wang <- mgeneSim(
  genes = arab_entrez,
  semData = sem_arab,
  measure = "Wang",
  combine = "BMA"
)

df_arab_wang <- as.data.frame(sim_arab_wang)

rownames(df_arab_wang) <- names(arab_entrez)
colnames(df_arab_wang) <- names(arab_entrez)

write.csv(df_arab_wang, file.path(output_dir, "wang_similarity_arabidopsis_genes.csv"))
cat("✔ Arabidopsis genes Wang matrix guardada con AGI\n")
####################################### PRUEBA 3 - Cálculo de similitud semántica entre genes de la arabdosis Thaliana.