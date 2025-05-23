# ClusteringSupervisado

Este paquete en R implementa dos algoritmos de *clustering supervisado* diseñados para seleccionar prototipos representativos por clase y clasificar nuevas observaciones en función de su proximidad a dichos prototipos.

###  Algoritmos incluidos

- **`supervised_pam()`**  
  Aplica el algoritmo PAM (Partitioning Around Medoids) de forma supervisada, seleccionando prototipos dentro de cada clase conocida.

- **`supervised_sridhcr()`**  
  Algoritmo iterativo inspirado en SRIDHCR que selecciona un prototipo aleatorio por clase y lo optimiza mediante búsqueda local (hill climbing) para maximizar la precisión de clasificación.

---

##  Instalación

Para instalar directamente desde GitHub:

```r
install.packages("devtools")  # si aún no lo tienes instalado
devtools::install_github("sarasanchezsanclemente/ClusteringSupervisado")
