ColorMixPlot <- function(object, 
                         features, 
                         slot = "data", 
                         colors = c("red", "green"), 
                         pt_size = 0.5,
                         section_number,
                         max_feature1 = NULL,
                         max_feature2 = NULL
) {
  stopifnot(length(features) == 2,
            is.character(features),
            length(colors) == 2,
            is.character(colors))
  
  # Subset object for sectionID of interest to make function faster
  object_sub <- object[, object$sectionID == section_number]
  
  # Retrieve data from the subset in the slot data of the current seurat object
  exprVals <- FetchData(object_sub, vars = features, slot = "data")
  
  # Apply max cutoff if specified in the command
  if (!is.null(max_feature1)) {
    exprVals[exprVals[, 1] > max_feature1, 1] <- max_feature1
  }
  if (!is.null(max_feature2)) {
    exprVals[exprVals[, 2] > max_feature2, 2] <- max_feature2
  }
  
  # Rescale the retrived values to get them between 0 and 1 in order to attribute color values
  exprVals_rescaled <- apply(exprVals, 2, scales::rescale)
  
  # Create colors for selected features
  ftr1_cols <- scales::gradient_n_pal(colours = c("#FFFFFF", colors[1], paste0("dark", colors[1])), values = c(0, 1))(exprVals_rescaled[, 1])
  ftr2_cols <- scales::gradient_n_pal(colours = c("#FFFFFF", colors[2], paste0("dark", colors[2])), values = c(0, 1))(exprVals_rescaled[, 2])
  
  # Create the color for overlapped points
  ftrs_blended_cols <- apply(cbind(ftr1_cols, ftr2_cols), 1, function(x) {
    # Check if either color is #FFFFFF (white)
    if (x[1] == "#FFFFFF") {
      blended_color <- x[2]
    } else if (x[2] == "#FFFFFF") {
      blended_color <- x[1]
    } else {
      # Blend colors when neither is #FFFFFF (white)
      blended_color <- colorjam::blend_colors(x, preset = "none")
    }
    return(blended_color)
  })
  
  # Create color list and use names of the features parameter
  cols_list <- list(ftr1_cols, ftr2_cols, ftrs_blended_cols) |> 
    setNames(nm = c(features, "mixed"))
  features <- features |> setNames(nm = features)
  
  # Create color gradients
  gradient_red <- c("#FFFFFF", colors[1], paste0("dark", colors[1]))
  gradient_green <- c("#FFFFFF", colors[2], paste0("dark", colors[2]))
  gradient_null <- "white"
  gradient_colors <- list(gradient_red, gradient_green, gradient_null) |> 
    setNames(nm = c(features, "mixed"))
  
  # Fetch the coordinate of each read to plot in 2D, subset to coordinate corresponding to the sectionID selected
  coor <- GetCoordinates(object)
  coor2 <- coor[coor$sampleID == section_number,]
  
  # Create ggplot object to plot
  gg <- bind_cols(exprVals, coor2)
  
  # Get maximum dimension from section_number to plot
  full_width <- object@tools[["Staffli"]]@image_info[["full_width"]][section_number]
  full_height <- object@tools[["Staffli"]]@image_info[["full_height"]][section_number]
  
  # Create list of the df exprVals, to efficiently use in the loop
  exprVals_list <- as.list(exprVals)
  
  ### Create the cell type proportion plots for each feature AND for the overlap 
  plots <- lapply(names(cols_list), function(nm) {
    cols <- cols_list[[nm]] #use the color corresponding to the feature used
    cell_proportion <- exprVals_list[[nm]] #use the cell_proportion corresponding to the feature used
    col_scale <- gradient_colors[[nm]]
    
    suppressWarnings({
      p1 <- ggplot(gg, aes(pxl_col_in_fullres, pxl_row_in_fullres)) +
        geom_point(color = cols, size = pt_size, aes(fill = cell_proportion)) +
        
        scale_fill_gradientn(colours = col_scale, #use color
                             breaks = c(round(min(cell_proportion, na.rm = TRUE), digits = 2),
                                        round((max(cell_proportion, na.rm = TRUE)*4/4), digits = 2),
                                        round((max(cell_proportion, na.rm = TRUE)*3/4), digits = 2), #add 4 ticks based on the respective cell_proportion
                                        round((max(cell_proportion, na.rm = TRUE)*2/4), digits = 2),
                                        round((max(cell_proportion, na.rm = TRUE)*1/4), digits = 2)),
                             limits = c(round(min(cell_proportion, na.rm = TRUE), digits = 2),
                                        round(max(cell_proportion, na.rm = TRUE), digits = 2))
        ) +
        
        scale_y_reverse(limits = c(full_height, 0)) +
        scale_x_continuous(limits = c(0, full_width)) +
        
        labs(x = NULL, y = NULL) +
        
        theme(panel.background = element_rect(fill = "grey",
                                              colour = "grey",
                                              size = 0,
                                              linetype = "solid"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = NULL,
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              aspect.ratio = 1) +
        
        ggtitle(label = ifelse(nm == "mixed", paste(features[1], features[2], sep = "+"), features[nm]))
    })
  })
  plot1 <- plots[[1]]
  plot2 <- plots[[2]]
  plot3 <- plots[[3]]
  
  plots <- plot1 + plot2 + plot3 +
    plot_layout(guides = "collect")
  
  return(plots)
}




ColorMixPlot3 <- function(object, 
                          features, 
                          slot = "data", 
                          colors = c("red", "green", "blue"), 
                          pt_size = 0.5,
                          section_number,
                          max_feature1 = NULL,
                          max_feature2 = NULL,
                          max_feature3 = NULL
) {
  stopifnot(length(features) == 3,
            is.character(features),
            length(colors) == 3,
            is.character(colors))
  
  # Subset object for sectionID of interest to make function faster
  object_sub <- object[, object$sectionID == section_number]
  
  # Retrieve data from the subset in the slot data of the current seurat object
  exprVals <- FetchData(object_sub, vars = features, slot = "data")
  
  # Apply max cutoff if specified in the command
  if (!is.null(max_feature1)) {
    exprVals[exprVals[, 1] > max_feature1, 1] <- max_feature1
  }
  if (!is.null(max_feature2)) {
    exprVals[exprVals[, 2] > max_feature2, 2] <- max_feature2
  }
  if (!is.null(max_feature3)) {
    exprVals[exprVals[, 3] > max_feature3, 3] <- max_feature3
  }
  
  # Rescale the retrived values to get them between 0 and 1 in order to attribute color values
  exprVals_rescaled <- apply(exprVals, 2, scales::rescale)
  
  # Create colors for selected features
  ftr1_cols <- scales::gradient_n_pal(colours = c("#FFFFFF", colors[1], paste0("dark", colors[1])), values = c(0, 1))(exprVals_rescaled[, 1])
  ftr2_cols <- scales::gradient_n_pal(colours = c("#FFFFFF", colors[2], paste0("dark", colors[2])), values = c(0, 1))(exprVals_rescaled[, 2])
  ftr3_cols <- scales::gradient_n_pal(colours = c("#FFFFFF", colors[3], paste0("dark", colors[3])), values = c(0, 1))(exprVals_rescaled[, 3])
  
  # Create the color for overlapped points
  ftrs_blended_cols <- apply(cbind(ftr1_cols, ftr2_cols, ftr3_cols), 1, function(x) {
    # Check if any color is #FFFFFF (white)
    if (x[1] == "#FFFFFF" && x[2] == "#FFFFFF") {
      blended_color <- x[3]
    } else if (x[1] == "#FFFFFF" && x[3] == "#FFFFFF") {
      blended_color <- x[2]
    } else if (x[2] == "#FFFFFF" && x[3] == "#FFFFFF") {
      blended_color <- x[1]
    } else {
      # Blend colors when none of them is #FFFFFF (white)
      blended_color <- colorjam::blend_colors(x, preset = "none")
    }
    return(blended_color)
  })
  
  # Create color list and use names of the features parameter
  cols_list <- list(ftr1_cols, ftr2_cols, ftr3_cols, ftrs_blended_cols) |> 
    setNames(nm = c(features, "mixed"))
  features <- features |> setNames(nm = features)
  
  # Create color gradients
  gradient_red <- c("#FFFFFF", colors[1], paste0("dark", colors[1]))
  gradient_green <- c("#FFFFFF", colors[2], paste0("dark", colors[2]))
  gradient_blue <- c("#FFFFFF", colors[3], paste0("dark", colors[3]))
  gradient_null <- "white"
  gradient_colors <- list(gradient_red, gradient_green, gradient_blue, gradient_null) |> 
    setNames(nm = c(features, "mixed"))
  
  # Fetch the coordinate of each read to plot in 2D, subset to coordinate corresponding to the sectionID selected
  coor <- GetCoordinates(object)
  coor2 <- coor[coor$sampleID == section_number,]
  
  # Create ggplot object to plot
  gg <- bind_cols(exprVals, coor2)
  
  # Get maximum dimension from section_number to plot
  full_width <- object@tools[["Staffli"]]@image_info[["full_width"]][section_number]
  full_height <- object@tools[["Staffli"]]@image_info[["full_height"]][section_number]
  
  # Create list of the df exprVals, to efficiently use in the loop
  exprVals_list <- as.list(exprVals)
  
  ### Create the cell type proportion plots for each feature AND for the overlap 
  plots <- lapply(names(cols_list), function(nm) {
    cols <- cols_list[[nm]] #use the color corresponding to the feature used
    cell_proportion <- exprVals_list[[nm]] #use the cell_proportion corresponding to the feature used
    col_scale <- gradient_colors[[nm]]
    
    suppressWarnings({
      p1 <- ggplot(gg, aes(pxl_col_in_fullres, pxl_row_in_fullres)) +
        geom_point(color = cols, size = pt_size, aes(fill = cell_proportion)) +
        
        scale_fill_gradientn(colours = col_scale, #use color
                             breaks = c(round(min(cell_proportion, na.rm = TRUE), digits = 2),
                                        round((max(cell_proportion, na.rm = TRUE)*4/4), digits = 2),
                                        round((max(cell_proportion, na.rm = TRUE)*3/4), digits = 2), #add 4 ticks based on the respective cell_proportion
                                        round((max(cell_proportion, na.rm = TRUE)*2/4), digits = 2),
                                        round((max(cell_proportion, na.rm = TRUE)*1/4), digits = 2)),
                             limits = c(round(min(cell_proportion, na.rm = TRUE), digits = 2),
                                        round(max(cell_proportion, na.rm = TRUE), digits = 2))
        ) +
        
        scale_y_reverse(limits = c(full_height, 0)) +
        scale_x_continuous(limits = c(0, full_width)) +
        
        labs(x = NULL, y = NULL) +
        
        theme(panel.background = element_rect(fill = "grey",
                                              colour = "grey",
                                              size = 0,
                                              linetype = "solid"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = NULL,
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              aspect.ratio = 1) +
        
        ggtitle(label = ifelse(nm == "mixed", paste(features[1], features[2], features[3], sep = "+"), features[nm]))
    })
  })
  plot1 <- plots[[1]]
  plot2 <- plots[[2]]
  plot3 <- plots[[3]]
  plot4 <- plots[[4]]
  
  plots <- plot1 + plot2 + plot3 + plot4 +
    plot_layout(guides = "collect", ncol = 4)
  
  return(plots)
}