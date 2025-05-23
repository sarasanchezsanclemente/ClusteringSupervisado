#' Algoritmo PAM supervisado
#'
#' Aplica el método PAM de forma supervisada seleccionando prototipos globales que minimizan la impureza y penalización.
#' según lo descrito en el artículo de Eick et al. (2004).
#' @param datos Matriz o data.frame de muestras.
#' @param clases Vector de clases.
#' @param k Número total de prototipos a seleccionar.
#' @param beta Peso de la penalización.
#' @param max_iter Número máximo de iteraciones para la fase de mejora.
#' @return Lista con prototipos, predicciones, impureza, coste y matriz de confusión.
#' @export
#' @import cluster
supervised_pam <- function(datos, clases, k = 2, beta = 1, max_iter = 50) {
  datos <- scale(datos)
  clases <- as.factor(clases)
  n <- nrow(datos)
  c <- length(levels(clases))

  penalty <- function(k, c, n) {
    if (k < c) return(9999)
    if (k == c) return(0)
    sqrt(k - c) / n
  }

  calcular_impureza <- function(asignaciones, clases) {
    impureza <- 0
    for (cl in unique(asignaciones)) {
      idx <- which(asignaciones == cl)
      if (length(idx) == 0) next
      mayoritaria <- names(which.max(table(clases[idx])))
      impureza <- impureza + sum(clases[idx] != mayoritaria)
    }
    impureza / length(clases)
  }

  clase_freq <- names(which.max(table(clases)))
  idx_inicial <- which(clases == clase_freq)[1]
  medoids_idx <- idx_inicial

  while (length(medoids_idx) < k) {
    candidatos <- setdiff(1:n, medoids_idx)
    mejor_q <- Inf
    mejor_candidato <- NULL
    for (cand in candidatos) {
      prueba_medoids <- c(medoids_idx, cand)
      asignaciones <- apply(datos, 1, function(x) {
        which.min(apply(datos[prueba_medoids, , drop = FALSE], 1, function(m) sum((x - m)^2)))
      })
      impureza <- calcular_impureza(asignaciones, clases)
      q <- impureza + beta * penalty(length(prueba_medoids), c, n)
      if (q < mejor_q) {
        mejor_q <- q
        mejor_candidato <- cand
      }
    }
    medoids_idx <- c(medoids_idx, mejor_candidato)
  }

  iter <- 0
  mejora <- TRUE
  while (mejora && iter < max_iter) {
    mejora <- FALSE
    iter <- iter + 1
    for (m in medoids_idx) {
      candidatos <- setdiff(1:n, medoids_idx)
      for (cand in candidatos) {
        nuevos_medoids <- medoids_idx
        nuevos_medoids[nuevos_medoids == m] <- cand
        asignaciones <- apply(datos, 1, function(x) {
          which.min(apply(datos[nuevos_medoids, , drop = FALSE], 1, function(m) sum((x - m)^2)))
        })
        impureza <- calcular_impureza(asignaciones, clases)
        q <- impureza + beta * penalty(length(nuevos_medoids), c, n)

        asign_base <- apply(datos, 1, function(x) {
          which.min(apply(datos[medoids_idx, , drop = FALSE], 1, function(m) sum((x - m)^2)))
        })
        impureza_base <- calcular_impureza(asign_base, clases)
        if (q < (impureza_base + beta * penalty(k, c, n))) {
          medoids_idx <- nuevos_medoids
          mejora <- TRUE
        }
      }
    }
  }

  prototipos <- datos[medoids_idx, , drop = FALSE]
  asignaciones <- apply(datos, 1, function(x) {
    which.min(apply(prototipos, 1, function(m) sum((x - m)^2)))
  })
  cl_to_class <- rep(NA, k)
  for (cl in 1:k) {
    idx <- which(asignaciones == cl)
    if (length(idx) > 0) {
      mayoritaria <- names(which.max(table(clases[idx])))
      cl_to_class[cl] <- mayoritaria
    }
  }

  predicciones <- cl_to_class[asignaciones]
  impureza_final <- calcular_impureza(asignaciones, clases)
  coste_total <- impureza_final + beta * penalty(k, c, n)
  precision <- mean(predicciones == clases)
  clases_predichas <- factor(predicciones, levels = levels(clases))

  return(list(
    medoids_idx = medoids_idx,
    medoids = prototipos,
    clustering = asignaciones,
    impureza = impureza_final,
    coste = coste_total,
    precision = precision,
    clases_predichas = clases_predichas,
    matriz_confusion = table(Real = clases, Prediccion = clases_predichas)
  ))
}


#' Algoritmo SRIDHCR supervisado
#'
#' Implementa el algoritmo de clustering supervisado SRIDHCR (Single Representative Insertion/Deletion Hill Climbing with Random Restart)
#' según lo descrito en el artículo de Eick et al. (2004).
#'
#' @param datos Matriz o data.frame de muestras (se estandariza internamente).
#' @param clases Vector de clases (factor).
#' @param beta Peso de la penalización por número de clusters.
#' @param max_iter Número máximo de iteraciones por reinicio.
#' @param reinicios Número de reinicios aleatorios para evitar mínimos locales.
#' @return Una lista con prototipos, predicciones, matriz de confusión, impureza, precisión y valor de la función q.
#' @export
supervised_sridhcr <- function(datos, clases, beta = 0.4, max_iter = 50, reinicios = 5) {
  datos <- scale(datos)
  clases <- as.factor(clases)
  n <- nrow(datos)
  c <- length(levels(clases))

  # Función de impureza: proporción de elementos mal clasificados en cada cluster
  calcular_impureza <- function(asignaciones, clases) {
    sum(sapply(unique(asignaciones), function(cl) {
      idx <- which(asignaciones == cl)
      if (length(idx) == 0) return(0)
      mayoritaria <- names(which.max(table(clases[idx])))
      sum(clases[idx] != mayoritaria)
    })) / length(clases)
  }

  # Penalización por número de clusters
  penalizacion <- function(k, c, n) {
    if (k < c) return(9999)
    if (k == c) return(0)
    sqrt(k - c) / n
  }

  # Función de coste q(X)
  evaluar_q <- function(indices_prototipos) {
    prototipos <- datos[indices_prototipos, , drop = FALSE]
    asignaciones <- apply(datos, 1, function(x) {
      which.min(apply(prototipos, 1, function(p) sum((x - p)^2)))
    })
    imp <- calcular_impureza(asignaciones, clases)
    q <- imp + beta * penalizacion(length(indices_prototipos), c, n)
    return(list(q = q, impureza = imp, asignaciones = asignaciones))
  }

  mejor_q <- Inf
  mejor_resultado <- NULL

  for (r in 1:reinicios) {
    k_inicial <- sample((c + 1):(2 * c), 1)
    actuales <- sample(1:n, k_inicial)

    iter <- 0
    continuar <- TRUE

    while (continuar && iter < max_iter) {
      iter <- iter + 1
      vecinos <- list()

      # Añadir un nuevo representante
      for (nuevo in setdiff(1:n, actuales)) {
        candidato <- c(actuales, nuevo)
        eval <- evaluar_q(candidato)
        vecinos[[length(vecinos) + 1]] <- list(indices = candidato, q = eval$q, imp = eval$impureza, asign = eval$asignaciones)
      }

      # Eliminar un representante
      if (length(actuales) > 1) {
        for (quitar in actuales) {
          candidato <- setdiff(actuales, quitar)
          eval <- evaluar_q(candidato)
          vecinos[[length(vecinos) + 1]] <- list(indices = candidato, q = eval$q, imp = eval$impureza, asign = eval$asignaciones)
        }
      }

      q_actual <- evaluar_q(actuales)$q
      mejor_vecino <- vecinos[[which.min(sapply(vecinos, function(v) v$q))]]

      if (mejor_vecino$q < q_actual || (mejor_vecino$q == q_actual && length(mejor_vecino$indices) > length(actuales))) {
        actuales <- mejor_vecino$indices
        asign_actual <- mejor_vecino$asign
        imp_actual <- mejor_vecino$imp
        q_actual <- mejor_vecino$q
      } else {
        continuar <- FALSE
      }
    }

    # Guardar la mejor solución global
    if (exists("q_actual") && q_actual < mejor_q) {
      prototipos <- datos[actuales, , drop = FALSE]
      predicciones <- apply(datos, 1, function(x) {
        dists <- apply(prototipos, 1, function(p) sum((x - p)^2))
        asign <- which.min(dists)
        mayoritaria <- names(which.max(table(clases[asign_actual == asign])))
        return(mayoritaria)
      })

      mejor_resultado <- list(
        prototipos = prototipos,
        indices_prototipos = actuales,
        clustering = asign_actual,
        clases_predichas = factor(predicciones, levels = levels(clases)),
        matriz_confusion = table(Real = clases, Prediccion = factor(predicciones, levels = levels(clases))),
        precision = mean(predicciones == clases),
        impureza = imp_actual,
        q = q_actual
      )
      mejor_q <- q_actual
    }
  }

  return(mejor_resultado)
}
