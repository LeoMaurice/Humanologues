library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(rrcov)
library(dbscan)
library(cluster)
library(ggrepel)
library(ggprism)

is.rrcovPca <- function(pca_res){
  return(
    class(pca_res) == "PcaHubert" | 
      class(pca_res) == "PcaCov" | 
      class(pca_res) == "PcaGrid" |
      class(pca_res) == "PcaLocantore" |
      class(pca_res) == "PcaProj"
  )
}

# Fonction pour les sorties standard
pca_visualization <- function(pca_res, docmetadatas, projected_data, log = FALSE) {
  
  # Helper functions to extract data based on PCA object class
  get_individual_coords <- function(pca_res) {
    if ("PCA" %in% class(pca_res)) {
      # FactoMineR PCA object
      ind_coords <- pca_res$ind$coord
    } else if (is.rrcovPca(pca_res)) {
      # rrcov PCA object
      ind_coords <- pca_res@scores
    } else {
      stop(paste("Unsupported PCA object of class", class(pca_res)))
    }
    return(ind_coords)
  }
  
  get_variable_coords <- function(pca_res) {
    if ("PCA" %in% class(pca_res)) {
      var_coords <- pca_res$var$coord
    } else if (is.rrcovPca(pca_res)) {
      var_coords <- pca_res@loadings
    } else {
      stop(paste("Unsupported PCA object of class", class(pca_res)))
    }
    return(var_coords)
  }
  
  get_eigenvalues <- function(pca_res) {
    if ("PCA" %in% class(pca_res)) {
      eig <- pca_res$eig
    } else if (is.rrcovPca(pca_res)) {
      eig <- pca_res@eigenvalues
    } else {
      stop(paste("Unsupported PCA object of class", class(pca_res)))
    }
    return(eig)
  }
  
  # Projection des individus sur l'espace PCA
  proj_individuals <- function(pca_res, docmetadatas, dims = c(1, 2)) {
    ind_coords <- get_individual_coords(pca_res)
    ind_data <- as.data.frame(ind_coords[, dims]) %>%
      rename_with(~ paste0("Dim", dims), everything()) 
    
    ind_data$document <- row.names(ind_data)
    
    ind_data <- ind_data%>%
      left_join(docmetadatas,
                by = "document")
    
    x_var <- paste0("PC", dims[1])
    y_var <- paste0("PC", dims[2])
    
    # Calculate percentage of variance explained
    eig <- get_eigenvalues(pca_res)
    if ("PCA" %in% class(pca_res)) {
      perc_var <- eig[, 2]
    } else if (is.rrcovPca(pca_res)) {
      perc_var <- eig / sum(eig) * 100
    }
    
    PCA_plot <- ggplot(projected_data, aes_string(x = x_var, y = y_var, color = "supp", label = "Nom")) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text_repel(
        box.padding = 0.5,  # Espace autour des labels
        max.overlaps = Inf, # Nombre maximal de chevauchements autorisés
        seed = 123          # Pour la reproductibilité
      ) +
      # geom_text(aes(label = Nom), vjust = -0.5, size = 5) +
      labs(
        title = paste("Projection sur les axes PCA", dims[1], "et", dims[2]),
        x = paste0("PCA ", dims[1], " (", round(perc_var[dims[1]], 1), "%)"),
        y = paste0("PCA ", dims[2], " (", round(perc_var[dims[2]], 1), "%)")
      ) +
      
      theme_pubclean() +
      theme(axis.text = element_text(size = 12),
            legend.text = element_text(20),
            axis.title = element_text(size = 20, face = "bold"),
            title = element_text(size = 26, face = "bold")) +
      scale_colour_prism(
        name = "Supplémentaire ?",
        palette = "colorblind_safe", 
        labels = c("Non", "Oui")
      )
    
    if (log) {
      PCA_plot <- PCA_plot +
        scale_x_continuous(trans = scales::pseudo_log_trans()) +
        scale_y_continuous(trans = scales::pseudo_log_trans()) +
        labs(title = paste("Projection sur les axes PCA", dims[1], "et", dims[2], "en échelle symlog"))
    }
    
    PCA_plot
  }
  
  # Pourcentage de représentation par axe
  plot_variance <- function(pca_res) {
    eig <- get_eigenvalues(pca_res)
    if ("PCA" %in% class(pca_res)) {
      eig_df <- as.data.frame(eig)
      eig_df$Axes <- seq_along(eig_df[, 1])
      eig_df$Variance <- eig_df[, 2]
    } else if (is.rrcovPca(pca_res)) {
      eig_df <- data.frame(
        eigenvalue = eig,
        Variance = eig / sum(eig) * 100,
        Cumulative = cumsum(eig) / sum(eig) * 100,
        Axes = seq_along(eig)
      )
    }
    
    ggplot(eig_df %>% slice(1:5), aes(x = Axes, y = Variance)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
      geom_text(aes(label = round(Variance, 1)), vjust = -0.5) +
      labs(
        title = "Pourcentage de variance expliquée par axe",
        x = "Axes",
        y = "Pourcentage de variance"
      ) +
      theme_pubclean()
  }
  
  # Variables qui comptent le plus sur chaque axe
  plot_top_variables <- function(pca_res, dim = 1, top_n = 25) {
    var_coords <- get_variable_coords(pca_res)
    
    var_data <- as.data.frame(var_coords[, dim, drop = FALSE]) %>%
      rename_with(~ paste0("Dim", dim), everything()) %>%
      mutate(Variable = rownames(var_coords)) %>%
      pivot_longer(cols = starts_with("Dim"), names_to = "Dimension", values_to = "Contribution") %>%
      group_by(Dimension) %>%
      slice_max(order_by = Contribution, n = top_n, with_ties = FALSE) %>%
      bind_rows(
        .,
        group_by(as.data.frame(var_coords[, dim, drop = FALSE]) %>%
                   rename_with(~ paste0("Dim", dim), everything()) %>%
                   mutate(Variable = rownames(var_coords)) %>%
                   pivot_longer(cols = starts_with("Dim"), names_to = "Dimension", values_to = "Contribution"),
                 Dimension) %>%
                   slice_min(order_by = Contribution, n = top_n, with_ties = FALSE)
      ) 
    
    ggplot(var_data, aes(x = reorder(Variable, Contribution), y = Contribution, fill = Dimension)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      labs(
        title = paste("Top", top_n, "variables contributrices (Dim", dim, ")"),
        x = "Variables",
        y = "Contribution"
      ) +
      theme_pubclean()
  }
  
  plot_cutoff_variables <- function(pca_res, dim = 1, cut_off = 0.1) {
    var_coords <- get_variable_coords(pca_res)
    
    var_data <- as.data.frame(var_coords[, dim, drop = FALSE]) %>%
      rename_with(~ paste0("Dim", dim), everything()) %>%
      mutate(Variable = rownames(var_coords)) %>%
      pivot_longer(cols = starts_with("Dim"), names_to = "Dimension", values_to = "Contribution") %>%
      group_by(Dimension) %>%
      filter(abs(Contribution) >= cut_off)
    
    ggplot(var_data, aes(x = reorder(Variable, Contribution), y = Contribution, fill = Dimension)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      labs(
        title = paste("variables contributrices (Dim", dim, ") avec une contribution absolue de plus de ", cut_off),
        x = "Variables",
        y = "Contribution"
      ) +
      theme_pubclean()
  }
  
  # Cercle des corrélations
  plot_correlation_circle <- function(pca_res, dims = c(1, 2)) {
    # fonctionnement pas testé, car pas utilisé
    var_coords <- get_variable_coords(pca_res)
    var_data <- as.data.frame(var_coords[, dims]) %>%
      rename_with(~ paste0("Dim", dims), everything()) %>%
      mutate(Variable = rownames(var_coords))
    
    # Calculate percentage of variance explained
    eig <- get_eigenvalues(pca_res)
    if ("PCA" %in% class(pca_res)) {
      perc_var <- eig[, 2]
    } else if (is.rrcovPca(pca_res)) {
      perc_var <- eig / sum(eig) * 100
    }
    
    # Unit circle data
    circle <- data.frame(
      x = cos(seq(0, 2 * pi, length.out = 100)),
      y = sin(seq(0, 2 * pi, length.out = 100))
    )
    
    x_var <- paste0("Dim", dims[1])
    y_var <- paste0("Dim", dims[2])
    
    ggplot() +
      geom_path(data = circle, aes(x = x, y = y), color = "gray") +
      geom_segment(
        data = var_data,
        aes_string(x = 0, y = 0, xend = x_var, yend = y_var),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "steelblue"
      ) +
      geom_text(
        data = var_data,
        aes_string(x = x_var, y = y_var, label = "Variable"),
        color = "darkred",
        hjust = 0.5, vjust = 1.5, size = 3
      ) +
      labs(
        title = paste("Cercle des corrélations (Axes", dims[1], "et", dims[2], ")"),
        x = paste0("PCA ", dims[1], " (", round(perc_var[dims[1]], 1), "%)"),
        y = paste0("PCA ", dims[2], " (", round(perc_var[dims[2]], 1), "%)")
      ) +
      xlim(-1.1, 1.1) + ylim(-1.1, 1.1) +
      coord_fixed() +
      theme_pubclean()
  }
  
  plot_cluster <-  function(pca_res, docmetadatas) {
    # Helper functions to extract data based on PCA object class
    get_individual_coords <- function(pca_res) {
      if ("PCA" %in% class(pca_res)) {
        # FactoMineR PCA object
        ind_coords <- pca_res$ind$coord
      } else if (is.rrcovPca(pca_res)) {
        # rrcov PCA object
        ind_coords <- pca_res@scores
      } else {
        stop(paste("Unsupported PCA object of class", class(pca_res)))
      }
      return(ind_coords)
    }
    
    get_variable_coords <- function(pca_res) {
      if ("PCA" %in% class(pca_res)) {
        var_coords <- pca_res$var$coord
      } else if (is.rrcovPca(pca_res)) {
        var_coords <- pca_res@loadings
      } else {
        stop(paste("Unsupported PCA object of class", class(pca_res)))
      }
      return(var_coords)
    }
    
    get_eigenvalues <- function(pca_res) {
      if ("PCA" %in% class(pca_res)) {
        eig <- pca_res$eig
      } else if (is.rrcovPca(pca_res)) {
        eig <- pca_res@eigenvalues
      } else {
        stop(paste("Unsupported PCA object of class", class(pca_res)))
      }
      return(eig)
    }
    
    ind_coords <- get_individual_coords(pca_res)
    
    dims = c(1,2,3)
    
    ind_data <- as.data.frame(ind_coords[, dims]) %>%
      rename_with(~ paste0("Dim", dims), everything()) 
    
    ind_data$document <- row.names(ind_data)
    
    ind_data <- ind_data%>%
      left_join(docmetadatas,
                by = "document")
    cah <- hclust(dist(ind_coords), method = "ward.D2")
    clusters <- cutree(cah, k = 4)
    
    ind_data$cluster <- as.factor(clusters)
    
    print(clusters)
    
    x_var <- "Dim1"
    y_var <- "Dim2"
    
    # Calculate percentage of variance explained
    eig <- get_eigenvalues(pca_res)
    if ("PCA" %in% class(pca_res)) {
      perc_var <- eig[, 2]
    } else if (is.rrcovPca(pca_res)) {
      perc_var <- eig / sum(eig) * 100
    }
    
    PCA_plot12 <- ggplot(ind_data, aes_string(x = x_var, y = y_var, color = "cluster")) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text(aes(label = Nom), vjust = -0.5, size = 5) +
      labs(
        title = paste("Projection sur les axes PCA", dims[1], "et", dims[2]),
        x = paste0("PCA ", dims[1], " (", round(perc_var[dims[1]], 1), "%)"),
        y = paste0("PCA ", dims[2], " (", round(perc_var[dims[2]], 1), "%)")
      ) +
      theme_pubclean() +
      theme(legend.position = "none",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 20, face = "bold"),
            title = element_text(size = 26, face = "bold")) 
    
    x_var <- "Dim2"
    y_var <- "Dim3"
    
    # Calculate percentage of variance explained
    eig <- get_eigenvalues(pca_res)
    if ("PCA" %in% class(pca_res)) {
      perc_var <- eig[, 2]
    } else if (is.rrcovPca(pca_res)) {
      perc_var <- eig / sum(eig) * 100
    }
    
    PCA_plot23 <- ggplot(ind_data, aes_string(x = x_var, y = y_var, color = "cluster")) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text(aes(label = Nom), vjust = -0.5, size = 5) +
      labs(
        title = paste("Projection sur les axes PCA", dims[1], "et", dims[2]),
        x = paste0("PCA ", dims[1], " (", round(perc_var[dims[1]], 1), "%)"),
        y = paste0("PCA ", dims[2], " (", round(perc_var[dims[2]], 1), "%)")
      ) +
      theme_pubclean() +
      theme(legend.position = "none",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 20, face = "bold"),
            title = element_text(size = 26, face = "bold")) 
    return(
      list(
        dims12 = PCA_plot12,
        dims12 = PCA_plot23
      )
    )
  }
  
  plots_cluster <- plot_cluster(pca_res, docmetadatas)
  
  # Génération des graphiques
  plots <- list(
    proj_1_2 = proj_individuals(pca_res, docmetadatas, dims = c(1, 2)),
    proj_1_3 = proj_individuals(pca_res, docmetadatas, dims = c(1, 3)),
    proj_2_3 = proj_individuals(pca_res, docmetadatas, dims = c(2, 3)),
    plot_circle_1_2 = plot_correlation_circle(pca_res, dims = c(1, 2)),
    plot_circle_1_3 = plot_correlation_circle(pca_res, dims = c(1, 3)),
    plot_circle_2_3 = plot_correlation_circle(pca_res, dims = c(2, 3)),
    variance = plot_variance(pca_res),
    top_vars_dim1 = plot_cutoff_variables(pca_res, dim = 1,cut_off = 0.1),
    top_vars_dim2 = plot_cutoff_variables(pca_res, dim = 2,cut_off = 0.1),
    top_vars_dim3 = plot_cutoff_variables(pca_res, dim = 3,cut_off = 0.1),
    cluster12 = plots_cluster$dims12,
    cluster23 = plots_cluster$dims23
  )
  
  return(plots)
}

get_contributions <- function(pca_res){
  if ("PCA" %in% class(pca_res)) {
    correlation_dimension <- as.data.frame(pca_res_humanologues$var$coord)
  } else if (is.rrcovPca(pca_res)) {
    correlation_dimension <- as.data.frame(pca_res@loadings)
  }
  correlation_dimension$word <- rownames(correlation_dimension)
  correlation_dimension
}
