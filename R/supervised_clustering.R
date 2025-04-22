#' Algoritme PAM supervisat
#'
#' Aplica el mètode PAM de forma supervisada seleccionant prototips per classe.
#'
#' @param datos Matriu o data.frame de mostres.
#' @param clases Vector de classes.
#' @param k_por_clase Nombre de prototips per classe.
#' @param lambda No utilitzat actualment.
#' @return Llista amb prototips, classes, prediccions, matriu de confusió i precisió.
#' @export
#' @import cluster
supervised_pam <- function(datos, clases, k_por_clase = 1, lambda = 0) {
  if (!is.data.frame(datos) && !is.matrix(datos)) {
    stop("'datos' ha de ser un data.frame o una matriu numèrica.")
  }
  if (length(clases) != nrow(datos)) {
    stop("La longitud de 'clases' ha de coincidir amb el nombre de files de 'datos'.")
  }
  clases <- as.factor(clases)
  datos <- scale(datos)
  medoides <- list()
  for (clase in levels(clases)) {
    datos_clase <- datos[clases == clase, , drop = FALSE]
    k <- min(k_por_clase, nrow(datos_clase))
    pam_result <- pam(datos_clase, k)
    indices_medoides <- pam_result$id.med
    medoides_clase <- datos_clase[indices_medoides, , drop = FALSE]
    medoides[[clase]] <- medoides_clase
  }
  prototipos <- do.call(rbind, medoides)
  prototipos_clases <- rep(names(medoides), each = k_por_clase)
  predecir_clase <- function(x) {
    distancias <- apply(prototipos, 1, function(p) sum((x - p)^2))
    clase_predicha <- prototipos_clases[which.min(distancias)]
    return(clase_predicha)
  }
  clases_predichas <- apply(datos, 1, predecir_clase)
  matriz_confusion <- table(Real = clases, Predicho = clases_predichas)
  precision <- mean(clases == clases_predichas)
  return(list(
    prototipos = prototipos,
    clases_prototipos = prototipos_clases,
    clases_predichas = clases_predichas,
    matriz_confusion = matriz_confusion,
    precision = precision,
    predecir = predecir_clase
  ))
}

#' Algoritme SRIDHCR supervisat
#'
#' Aplica l'algoritme SRIDHCR mitjançant optimització iterativa.
#'
#' @param datos Matriu o data.frame de mostres.
#' @param clases Vector de classes.
#' @param max_iter Nombre màxim d'iteracions per reinici.
#' @param max_reinicios Nombre màxim de reinicis aleatoris.
#' @return Millor resultat obtingut (prototips, prediccions, precisió, etc.).
#' @export
supervised_sridhcr <- function(datos, clases, max_iter = 20, max_reinicios = 5) {
  clases <- as.factor(clases)
  datos <- scale(datos)
  muestras_por_clase <- split(1:nrow(datos), clases)
  mejor_precision <- -Inf
  mejor_resultado <- NULL
  for (reinicio in 1:max_reinicios) {
    prototipos <- sapply(muestras_por_clase, function(indices) sample(indices, 1))
    sin_mejora <- 0
    for (iter in 1:max_iter) {
      predichas <- apply(datos, 1, function(x) {
        distancias <- apply(datos[prototipos, , drop = FALSE], 1, function(p) sum((x - p)^2))
        names(prototipos)[which.min(distancias)]
      })
      precision_actual <- mean(predichas == clases)
      if (precision_actual > mejor_precision) {
        mejor_precision <- precision_actual
        mejor_resultado <- list(
          prototipos = datos[prototipos, , drop = FALSE],
          clases_prototipos = names(prototipos),
          clases_predichas = predichas,
          matriz_confusion = table(Real = clases, Predicho = predichas),
          precision = precision_actual,
          predecir = function(x) {
            distancias <- apply(mejor_resultado$prototipos, 1, function(p) sum((x - p)^2))
            mejor_resultado$clases_prototipos[which.min(distancias)]
          }
        )
        sin_mejora <- 0
      } else {
        sin_mejora <- sin_mejora + 1
      }
      for (clase in names(prototipos)) {
        candidatos <- muestras_por_clase[[clase]]
        for (nuevo_prot in candidatos) {
          if (nuevo_prot == prototipos[clase]) next
          nuevos_prot <- prototipos
          nuevos_prot[clase] <- nuevo_prot
          nuevas_pred <- apply(datos, 1, function(x) {
            distancias <- apply(datos[nuevos_prot, , drop = FALSE], 1, function(p) sum((x - p)^2))
            names(nuevos_prot)[which.min(distancias)]
          })
          nueva_precision <- mean(nuevas_pred == clases)
          if (nueva_precision > precision_actual) {
            prototipos <- nuevos_prot
            break
          }
        }
      }
      if (sin_mejora >= 5) break
    }
  }
  return(mejor_resultado)
}
