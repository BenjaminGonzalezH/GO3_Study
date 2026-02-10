# GO3 Study

Proyecto de estudio sobre Gene Ontology (GO) y medidas de similitud semántica. El repositorio contiene notebooks de revisión, un script en R para calcular matrices de similitud (Wang), datasets GO/GAF y resultados generados.

**Estado**: investigación/estudio personal, con insumos y resultados reproducibles.

## Contenido
- `1. Notas de revisión GO3.ipynb`: notas y revisión teórica.
- `2. Archivo a archivo (RustCode).ipynb`: exploración/implementación paso a paso.
- `3. Preguntas.ipynb`: preguntas y ejercicios de estudio.
- `R Estudio de Genes.R`: script principal en R para cálculos de similitud.
- `Files/`: archivos GO y GAF de referencia (por ejemplo, `go.obo`, `goa_human.gaf`, `tair.gaf`).
- `Graphs/`: gráficos exportados de benchmarks.
- `Images/`: imágenes de apoyo.
- `Papers/`: artículos base del estudio.
- `R_results/`: resultados CSV generados por el script en R.

## Requisitos
- R 4.x
- Paquetes Bioconductor: `GOSemSim`, `GO.db`, `org.Hs.eg.db`, `org.At.tair.db`

El script instala los paquetes si no están disponibles.

## Uso rápido (R)
Ejecuta el script para generar las matrices de similitud y guardar CSV en `R_results/`.

```r
# En R
source("R Estudio de Genes.R")
```

Salida esperada:
- `R_results/wang_similarity_GO_terms.csv`
- `R_results/wang_similarity_human_genes.csv`
- `R_results/wang_similarity_arabidopsis_genes.csv`

## Notas
- Los notebooks pueden requerir rutas locales o datos presentes en `Files/`.
- Los datasets GO/GAF son grandes; evita duplicarlos innecesariamente.
- Si cambias la ontología o el conjunto de genes en el script, los resultados se regeneran en `R_results/`.

## Licencia
This project is licensed under the MIT License – see the LICENSE file for details.
