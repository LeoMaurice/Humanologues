library(tidyverse)
library(FactoMineR)
library(factoextra)
library(rrcov)

# PCA data
prepepare_pca_data <- function(df){
  # Étape 1 : Mise en format large
  df_wide <- prepare_pca_data(df)
  
  # Étape 2 : Normalisation
  df_normalized <- normalize_frequencies(df_wide)
  
  return(df_normalized)
}

# Pipeline complet
perform_pca_pipeline <- function(df, df_normalized = data.frame(), robust = F, supplementaires) {
  # Etape 1 et 2
  if(is_empty(df_normalized)){
    df_normalized <- prepepare_pca_data(df)
  }
  
  # ids_supp <- df%>%
  #   filter(supp)%>%
  #   distinct(document)%>%
  #   pull(document)
  
  # Étape 3 : PCA
  pca_res <- perform_pca(df_normalized%>%
                           filter(!Nom %in% supplementaires), robust = robust)
  
  # Etape 4 : Projection
  df_supp_norm <- df_normalized%>%
    filter(Nom %in% supplementaires)
  new_data_centered <- sweep(df_supp_norm%>%
                               select(-Nom), 2, pca_res$center, "-")
  projected_data_supp <- as.matrix(new_data_centered) %*% as.matrix(pca_res@loadings)
  
  df_main_norm <- df_normalized%>%
    filter(!Nom %in% supplementaires)
  new_data_centered <- sweep(df_main_norm%>%
                               select(-Nom), 2, pca_res$center, "-")
  projected_data_main <- as.matrix(new_data_centered) %*% as.matrix(pca_res@loadings)
  
  projected_data <- bind_rows(
    as.data.frame(projected_data_main),
    as.data.frame(projected_data_supp)
  )%>%
    # mutate(doc_id = as.character(str_pad(row_number(), width = 2, pad = "0")))
    rownames_to_column("Nom")
  
  
  
  return(list(
    pca_res = pca_res,
    projected_data = projected_data
  ))
}
