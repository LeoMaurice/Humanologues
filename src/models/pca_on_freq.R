library(tidyverse)
library(FactoMineR)
library(factoextra)
library(Factoshiny)

# Étape 1 : Pivot Wider
prepare_pca_data <- function(df) {
  df_wide <- df %>%
    select(Nom,term,count)%>%
    pivot_wider(names_from = term, values_from = count, values_fill = 0, names_prefix = "freq_") %>%
    as.data.frame() # Mettre à 0 les valeurs manquantes
  
  rownames(df_wide) <- df_wide$Nom
  return(df_wide)
}

# Étape 2 : Normalisation des fréquences
normalize_frequencies <- function(df_wide) {
  # Exclure la colonne "document" des normalisations
  df_wide <- df_wide %>%
    mutate(across(starts_with("freq_"), scale))  # Normalisation (z-score) sur les colonnes "freq_"
  return(df_wide)
}

# Étape 3 : Réalisation de la PCA
perform_pca <- function(df_wide, robust = F) {
  # Exclure la colonne "document" pour la PCA
  if(robust)
    pca_res <- PcaHubert(df_wide%>%select(starts_with("freq_")), k = 5)
  else
    pca_res <- PCA(df_wide%>%select(starts_with("freq_")), graph = FALSE)
  return(pca_res)
}

# Étape 4 : Affichage interactif
display_pca_interactive <- function(pca_res, df_wide) {
  # Associer les noms des documents pour l'interactivité
  # pca_res$ind$coord <- cbind(documents = df_wide$Nom, pca_res$ind$coord)
  
  # Lancer Factoshiny pour un affichage interactif
  Factoshiny(pca_res)
}