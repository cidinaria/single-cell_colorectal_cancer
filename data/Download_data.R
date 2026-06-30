# Caso não tenha o BiocManager instalado
# install.packages("BiocManager")
# BiocManager::install("GEOquery")

library(GEOquery)

# Baixar o objeto do GEO
gse <- getGEO("GSE236581", GSEMatrix = TRUE)

# Se houver mais de uma plataforma, gse será uma lista. 
# Para pegar o primeiro elemento (GSEExpressionSet):
eSet <- gse[[1]]

# Obter a matriz de expressão
expressao <- exprs(eSet)

# Obter os metadados (dados clínicos/fenotípicos dos pacientes)
metadados <- pData(eSet)
