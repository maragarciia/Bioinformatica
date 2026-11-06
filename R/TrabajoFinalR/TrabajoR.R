##########################################################################################
#####################                                        		##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       		##################
##########################################################################################

Mara García García, NP: 153711

### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

datos <- read.table("archivo-entrada.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

head(datos)
tail(datos)
dim(datos)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

control <- datos[, c("C1-TPM", "C2-TPM", "C3-TPM")]
condicion <- datos[, c("1W5-TPM", "1W6-TPM", "1W5-COV")]

datos$Media_Control <- apply(control, 1, mean)
datos$DE_Control <- apply(control, 1, sd)
datos$Media_Condicion <- apply(condicion, 1, mean)
datos$DE_Condicion <- apply(condicion, 1, sd)

head(datos)

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

col_ctrl <- "skyblue"
col_cond <- "orange"

par(mfrow = c(1, 2))
hist(datos$Media_Control,   col = col_ctrl, main = "Media Control",
     xlab = "Expresión media (TPM)", ylab = "Frecuencia")
hist(datos$Media_Condicion, col = col_cond, main = "Media Condición",
     xlab = "Expresión media (TPM)", ylab = "Frecuencia")
par(mfrow = c(1, 1))

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

boxplot(datos$Media_Control, datos$Media_Condicion,
        names = c("Control", "Condición"),
        col = c(col_ctrl, col_cond),
        main = "Expresión: Control vs Condición",
        ylab = "Expresión media (TPM)")
legend("topleft", legend = c("Control", "Condición"), fill = c(col_ctrl, col_cond))

### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

shapiro.test(datos$Media_Control)
shapiro.test(datos$Media_Condicion)

### Al evaluar la normalidad mediante la prueba de Shapiro-Wilk, se descarta la distribución normal 
### en los datos. Debido a esto, se aplica la prueba no paramétrica de Wilcoxon gen a gen.

p_values <- numeric(nrow(datos))

for (i in 1:nrow(datos)) {
  x <- as.numeric(control[i, ])
  y <- as.numeric(condicion[i, ])
  p_values[i] <- suppressWarnings(wilcox.test(x, y)$p.value)
}

datos$pvalue <- p_values

summary(datos$pvalue)   # no debe haber NA
min(datos$pvalue)       # p-value mínimo
sum(datos$pvalue < 0.05)    # debe dar 0

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

### Con un tamaño muestral de n=3 por grupo, el p-valor mínimo alcanzable estadísticamente por 
### esta prueba es 0.1, por lo que no se obtienen genes con p < 0.05.
### Por ello, se sube el umbral a pvalue <= 0.1.

deg <- datos[datos$pvalue <= 0.1, ]

genes_deg <- rownames(deg)
total_deg <- length(genes_deg)

sobreexpresados <- rownames(deg[deg$Media_Condicion > deg$Media_Control, ])
reprimidos      <- rownames(deg[deg$Media_Condicion < deg$Media_Control, ])

cat("Número total de genes diferencialmente expresados (p <= 0.1):", total_deg, "\n")
cat("Genes sobreexpresados (Condición > Control):", length(sobreexpresados), "\n")
cat("Genes reprimidos (Condición < Control):", length(reprimidos), "\n")

### Resultado: 429 genes diferencialmente expresados (p <= 0.1), de los cuales 69 están
### sobreexpresados y 360 reprimidos en la Condición respecto al Control.

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

if (total_deg > 0) {
  # Extraer la matriz con las columnas de Control y Condición
  mat <- as.matrix(deg[, c(colnames(control), colnames(condicion))])
  
  heatmap(mat, scale = "row", Colv = NA, labRow = NA,
          col = colorRampPalette(c(col_ctrl, "white", col_cond))(50),
          ColSideColors = c(rep(col_ctrl, 3), rep(col_cond, 3)),
          main = "Heatmap de Genes DE (p <= 0.1)")
} else {
  cat("No se encontraron genes significativos para representar en el heatmap.\n")
}

### La barra superior indica el grupo (celeste = Control, naranja = Condición).
### En el mapa de calor, naranja = expresión alta y celeste = expresión baja de cada gen
### (escalado por filas).