server <- function(input, output, session) {
  
  # Wait 500 milliseconds after the user stops clicking before updating plots
  genes_debounced <- reactive({ input$abundance_genes }) %>% debounce(500)
  tools_debounced <- reactive({ input$tool_abundance }) %>% debounce(500)
  
  # Helper to translate the UI threshold (e.g., "0.0" or "60.0") to the list name ("default" or "60")
  get_thr_name <- function(thr_val) {
    thr_str <- as.character(thr_val)
    thr_str <- gsub("\\.0$", "", thr_str) # This strips the trailing .0
    
    if (thr_str %in% c("0", "default")) {
      return("default")
    }
    return(thr_str)
  }
  
  ## Caching Unigenes for a fast loading
  unigenes_filtered <- reactive({
    req(input$threshold_unigenes_id)
    thr_str <- gsub("\\.0$", "", as.character(input$threshold_unigenes_id))
    thr <- as.numeric(thr_str)
    
    if (thr > 0) {
      data_list$unigenes_base %>% 
        dplyr::filter(!(tool %in% c("DeepARG", "RGI-DIAMOND") & id < thr))
    } else {
      data_list$unigenes_base
    }
  })
  
  output$plot_count_genes_tool <- renderPlot({
    plot_count_genes_tool(
      unigenes = unigenes_filtered(), 
      tools_for_figure = input$tools_unigenes, 
      general_size = general_size, 
      pal_10_q = pal_10_q, 
      tool_label = tools_labels , 
      tools_levels = tools_levels, 
      texture = tools_texture) + 
      theme(panel.background = element_rect(colour = "black", fill = NA)) + 
      theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"))
  })
  
  output$plot_alluvial_classes <- renderPlot({
    plot_alluvial_classes(unigenes = unigenes_filtered(), 
                          levels_unigenes = data_list$levels_unigenes, 
                          threshold_plot = 0.99, remove_class_threshold = 0.005, 
                          tools_to_plot = input$tools_unigenes, 
                          tools_labels = tools_labels, 
                          tools_factors = tools_levels, 
                          pal_10_q = pal_10_q, 
                          general_size = general_size, 
                          gene_classes_list = gene_classes_list
    ) +
      theme(panel.background = element_rect(colour = "black", fill = NA))
  })
  
  # Pan- & Core- Resistome Plot
  
  pan_core <- reactive({
    req(input$threshold_pan_core_id, input$threshold_samples, input$threshold_proportion,
        input$tool_pan_core, input$environment_pan_core)
    
    thr_name <- get_thr_name(input$threshold_pan_core_id)
    
    sumpan2_df <- sumpan2_prepped[[thr_name]]
    core_df    <- core_prepped[[thr_name]]
    
    req(!is.null(sumpan2_df), !is.null(core_df))
    
    core_sum <- sum_core_adjust(
      core_df,
      as.numeric(input$threshold_samples),
      input$threshold_proportion
    ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(tool, habitat) %>%
      dplyr::summarise(core = sum(unigenes, na.rm = TRUE), .groups = "drop")
    
    sumpan2_df %>%
      dplyr::left_join(core_sum, by = c("tool", "habitat")) %>%
      dplyr::mutate(
        core = tidyr::replace_na(core, 0),
        prop = core / md,
        texture = ifelse(tool %in% tools_texture, "yes", "no")
      ) %>%
      dplyr::filter(
        tool %in% input$tool_pan_core,
        habitat %in% input$environment_pan_core
      ) %>%
      dplyr::mutate(
        tool = factor(as.character(tool),
                      levels = tools_levels[tools_levels %in% input$tool_pan_core])
      )
  }) %>% bindCache(
    input$threshold_pan_core_id,
    input$threshold_samples,
    input$threshold_proportion,
    sort(input$tool_pan_core),
    sort(input$environment_pan_core)
  )
  
  
  output$plot_pan_core_resistome <- renderPlot({
    req(input$tool_pan_core, input$environment_pan_core)
    df <- pan_core()
    req(nrow(df) > 0)
    
    tools_order <- match(input$tool_pan_core, tools_levels)
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]
    
    pal_figure <- pal_10_q[tools_order]
    tools_labels_figure <- tools_labels[tools_order]
    
    df %>%
      dplyr::select(!c(md, sd)) %>%
      tidyr::pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>%
      dplyr::mutate(
        metric = dplyr::recode(metric, mn = "Pan-resistome", core = "Core-resistome"),
        metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))
      ) %>%
      ggplot2::ggplot(ggplot2::aes(x = habitat, y = value)) +
      ggplot2::geom_jitter(
        ggplot2::aes(fill = tool, shape = texture),
        color = "black", stroke = 0.3, size = 2.5,
        width = 0.5, height = 0
      ) +
      ggplot2::facet_grid(metric ~ habitat, scales = "free") +
      ggplot2::scale_fill_manual(values = pal_figure, labels = lab_fn(tools_labels_figure), name = NULL) +
      ggplot2::scale_shape_manual(values = c(21, 24)) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          override.aes = list(shape = shape_tools, fill = pal_figure)
        ),
        shape = "none"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("") +
      ggplot2::ylab("ARGs") +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text.x   = ggplot2::element_text(size = general_size),
        legend.text    = ggplot2::element_text(size = general_size),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(15, 15, 15, 15, unit = "pt"),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = grid::unit(0, "pt"),
        title = ggplot2::element_text(size = general_size + 2, face = "bold"),
        axis.title = ggplot2::element_text(size = general_size + 1, face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = general_size),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(colour = "black", fill = NA)
      )
  })
  
  
  # Abundance Plot
  
  abundance_tool_sample_reactive <- reactive({
    req(input$threshold_abundance_id)
    thr_name <- get_thr_name(input$threshold_abundance_id)
    df <- abundance_prepped[[thr_name]]
    req(!is.null(df))
    df
  }) %>% bindCache(input$threshold_abundance_id)
  
  abundance_filtered <- reactive({
    req(tools_debounced(), input$environment_abundance)
    
    df <- abundance_tool_sample_reactive()
    
    df %>%
      dplyr::filter(
        tool %in% tools_debounced(),
        habitat %in% input$environment_abundance
      )
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance)
  )
  
  plot_abundance_reactive <- reactive({
    req(tools_debounced(), input$environment_abundance, input$threshold_abundance_id)
    
    # Lock the order
    sel_tools <- intersect(tools_levels, tools_debounced())
    
    df <- abundance_filtered() %>%
      dplyr::mutate(tool = factor(as.character(tool), levels = sel_tools))
    req(nrow(df) > 0)
    
    plot_total_abundance_diversity_new_version(
      dataset = df,
      tools_labels = tools_labels,
      tools_to_plot = tools_debounced(),       
      environments_plot = input$environment_abundance,
      general_size = general_size,
      pal_10_q = pal_10_q,
      metric = "abundance",
      sd = 2025,
      obs = 200,
      texture = tools_texture,
      tools_levels = tools_levels
    ) +
      theme(legend.position = "none")
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),                      
    sort(input$environment_abundance)
  )
  
  output$plot_abundance <- renderPlot(plot_abundance_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance)
  )
  output$plot_abundance_overview <- renderPlot(plot_abundance_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance)
  )
  
  plot_diversity_reactive <- reactive({
    req(tools_debounced(), input$environment_abundance, input$threshold_abundance_id)
    
    # Locking the order
    sel_tools <- intersect(tools_levels, tools_debounced())
    
    df <- abundance_filtered() %>%
      dplyr::mutate(tool = factor(as.character(tool), levels = sel_tools))
    req(nrow(df) > 0)
    
    plot_total_abundance_diversity_new_version(
      dataset = df,
      tools_labels = tools_labels,
      tools_to_plot = tools_debounced(),         
      environments_plot = input$environment_abundance,
      general_size = general_size,
      pal_10_q = pal_10_q,
      metric = "diversity",
      sd = 2025,
      obs = 200,
      texture = tools_texture,
      tools_levels = tools_levels
    ) +
      theme(legend.position = "none")
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),                      
    sort(input$environment_abundance)
  )
  
  output$plot_diversity <- renderPlot(plot_diversity_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance)
  )
  output$plot_diversity_overview <- renderPlot(plot_diversity_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance)
  )
  
  outputOptions(output, "plot_abundance", suspendWhenHidden = TRUE)
  outputOptions(output, "plot_abundance_overview", suspendWhenHidden = TRUE)
  outputOptions(output, "plot_diversity", suspendWhenHidden = TRUE)
  outputOptions(output, "plot_diversity_overview", suspendWhenHidden = TRUE)
  
  
  # Median abundance and diversity Plot
  
  abundance_class_reactice <- reactive({
    req(input$threshold_abundance_id, tools_debounced())
    thr_name <- get_thr_name(input$threshold_abundance_id)
    
    df <- abundance_class_prepped[[thr_name]]
    req(!is.null(df), nrow(df) > 0)
    
    sel_tools <- intersect(tools_levels, tools_debounced())
    
    df %>%
      dplyr::filter(tool %in% sel_tools) %>%                        
      dplyr::mutate(tool = factor(as.character(tool), levels = sel_tools))       
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced())
  )
  
  plot_abundance_class_reactive <- reactive({
    req(tools_debounced(), input$environment_abundance)
    
    df <- abundance_class_reactice()
    req(nrow(df) > 0)
    
    sel_tools <- intersect(tools_levels, tools_debounced())
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)
    
    plot_abundance_class_more_environments(
      df,
      input$environment_abundance, general_size,
      pal_sel,                      
      genes_debounced(),
      data_type = "abundance",
      other = input$plot_other,
      tools_levels = sel_tools,
      sel_tools,
      texture = tools_texture,
      pattern_density = pattern_density,
      pattern_spacing = pattern_spacing,
      pattern_fill = pattern_fill,
      pattern_size = pattern_size
    ) +
      theme(legend.position = "none")
  })
  
  output$plot_abundance_class <- renderPlot(plot_abundance_class_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance),
    sort(genes_debounced()),
    input$plot_other
  )
  output$plot_abundance_class_overview <- renderPlot(plot_abundance_class_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance),
    sort(genes_debounced()),
    input$plot_other
  )
  
  plot_diversity_class_reactive <- reactive({
    req(tools_debounced(), input$environment_abundance, genes_debounced())
    
    df <- abundance_class_reactice()
    req(nrow(df) > 0)
    
    sel_tools <- intersect(tools_levels, tools_debounced())
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)
    
    plot_abundance_class_more_environments(
      df,
      input$environment_abundance, general_size,
      pal_sel,                        
      genes_debounced(),
      data_type = "diversity",
      other = input$plot_other,
      tools_levels = sel_tools,
      sel_tools,
      texture = tools_texture,
      pattern_density = pattern_density,
      pattern_spacing = pattern_spacing,
      pattern_fill = pattern_fill,
      pattern_size = pattern_size
    ) +
      theme(legend.position = "none")
  })
  
  output$plot_diversity_class <- renderPlot(plot_diversity_class_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance),
    sort(genes_debounced()),
    input$plot_other
  )
  output$plot_diversity_class_overview <- renderPlot(plot_diversity_class_reactive()) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance),
    sort(genes_debounced()),
    input$plot_other
  )
  
  plot_diversity_class_legend <- reactive({
    req(input$threshold_abundance_id, input$environment_abundance, tools_debounced(), genes_debounced())
    
    df <- abundance_class_reactice()
    req(nrow(df) > 0)
    
    sel_tools <- intersect(tools_levels, tools_debounced())
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)
    
    p <- plot_abundance_class_more_environments(
      df,
      input$environment_abundance, general_size,
      pal_sel,                     
      genes_debounced(),
      data_type = "diversity",
      other = input$plot_other,
      tools_levels = sel_tools,
      sel_tools,
      texture = tools_texture,
      pattern_density = 0.01,
      pattern_spacing = 0.005,
      pattern_fill = "white",
      pattern_size = 0.4
    ) +
      theme(legend.position = "bottom")
    
    g_legend(p)
  })
  
  output$plot_diversity_class_legend <- renderPlot({
    grid::grid.newpage()
    grid::grid.draw(plot_diversity_class_legend())
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance),
    sort(genes_debounced()),
    input$plot_other
  )
  
  output$plot_diversity_class_legend_overview <- renderPlot({
    grid::grid.newpage()
    grid::grid.draw(plot_diversity_class_legend())
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance),
    sort(genes_debounced()),
    input$plot_other
  )
  
  plot_abundance_class_legend <- reactive({
    req(input$threshold_abundance_id, input$environment_abundance, tools_debounced(), genes_debounced())
    
    df <- abundance_class_reactice()
    req(nrow(df) > 0)
    
    sel_tools <- intersect(tools_levels, tools_debounced())
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)
    
    p <- plot_abundance_class_more_environments(
      df,
      input$environment_abundance,
      general_size,
      pal_sel,
      genes_debounced(),
      data_type = "abundance",
      other = input$plot_other,
      tools_levels = sel_tools,
      sel_tools,
      texture = tools_texture,
      pattern_density = 0.01,
      pattern_spacing = 0.005,
      pattern_fill = "white",
      pattern_size = 0.4
    ) +
      theme(legend.position = "bottom")
    
    g_legend(p)
  })
  
  output$plot_abundance_class_legend <- renderPlot({
    grid::grid.newpage()
    grid::grid.draw(plot_abundance_class_legend())
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(tools_debounced()),
    sort(input$environment_abundance),
    sort(genes_debounced()),
    input$plot_other
  )
  
  # Overlaps Plot
  
  recall_base <- reactive({
    req(input$threshold_overlap_id)
    thr <- get_thr_name(input$threshold_overlap_id)
    
    df <- switch(
      thr,
      "60" = data_list$recall_fnr60,
      "70" = data_list$recall_fnr70,
      "80" = data_list$recall_fnr80,
      data_list$recall_fnr
    )
    req(!is.null(df))
    
    df %>%
      dplyr::filter(
        new_level %in% input$overlap_genes,
        tool_ref %in% input$tool_overlap,
        tool_comp %in% input$tool_overlap_calc
      ) %>%
      dplyr::mutate(
        new_level = factor(new_level, levels = input$overlap_genes),
        tool_ref  = factor(tool_ref, levels = tools_levels),
        texture   = ifelse(tool_comp %in% tools_texture, "yes", "no")
      )
  }) %>%
    bindCache(
      input$threshold_overlap_id,
      input$overlap_genes,
      input$tool_overlap,
      input$tool_overlap_calc
    )
  
  
  plot_overlap_cstc  <- reactive({
    df <- recall_base()
    req(nrow(df) > 0)
    
    df %>%
      ggplot(aes(x = new_level, fill = tool_ref, y = recall)) +
      geom_violin() +
      ggrastr::rasterise(
        geom_jitter(aes(color = tool_comp, shape = texture), stroke = 1, size = 2.5, width = 0.1,  height = 0),
        dpi = 150
      ) + 
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                        labels = tools_levels[tools_levels %in% input$tool_overlap]) +
      scale_color_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
                         labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
      scale_x_discrete(labels = function(x) {
        x <- gsub("-", "-\n", x)
        x <- gsub(" ", "\n", x)
        x}) + 
      scale_shape_manual(values = c("no" = 16, "yes" = 17)) +
      facet_grid(tool_ref ~ new_level, scales = "free_x") +
      scale_y_continuous(limits = c(-0.1,1.1), 
                         breaks = seq(0, 1, length.out = 3),
                         labels = scales::label_number()) + 
      ylab("CSTC") + 
      xlab("") + 
      theme(
        legend.position = "none",
        legend.text = element_text(size = general_size ),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(15, 15, 15, 15, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        panel.background = element_rect(colour = "black", fill = NA),
        strip.text = element_blank(),
        strip.background = element_blank())
  })
  
  output$overlap_cstc <- renderPlot({plot_overlap_cstc()})
  output$overlap_cstc_overview <- renderPlot({plot_overlap_cstc()})
  
  plot_overlap_cstc_summary  <- reactive({
    df <- recall_base()
    req(nrow(df) > 0)
    
    df %>%
      filter(!is.na(recall)) %>% 
      ungroup() %>% 
      group_by(tool_ref, new_level) %>% 
      summarise(recall = median(recall), .groups = "drop") %>%
      ggplot(aes(x = "Class medians", fill = tool_ref, y = recall)) +
      geom_violin() +
      ggrastr::rasterise(
        geom_jitter(size = 1.5, width = 0.4, height = 0),
        dpi = 150
      ) + 
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                        labels = tools_levels[tools_levels %in% input$tool_overlap]) +
      facet_grid(tool_ref ~ ., scales = "free_x") +
      scale_y_continuous(limits = c(-0.2,1.2), 
                         breaks = seq(0, 1, length.out = 3),
                         labels = scales::label_number()) +
      scale_x_discrete(labels = function(x) {
        x <- gsub("-", "-\n", x)
        x <- gsub(" ", "\n", x)
        x}) + 
      ylab("") + 
      xlab("") + 
      theme(
        legend.position = "none",
        legend.text = element_text(size = general_size ),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(15, 15, 15, 15, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA),
        strip.text = element_text(size = general_size, face = "bold"),
        strip.background = element_blank())
  })
  
  output$overlap_cstc_summary <- renderPlot({plot_overlap_cstc_summary()})
  output$overlap_cstc_summary_overview <- renderPlot({plot_overlap_cstc_summary()})
  
  plot_overlap_csno  <- reactive({
    df <- recall_base()
    req(nrow(df) > 0)
    
    df %>%
      ggplot(aes(x = new_level, fill = tool_ref, y = fnr)) +
      geom_violin() +
      ggrastr::rasterise(
        geom_jitter(aes(color = tool_comp, shape = texture), stroke = 1, size = 2.5, width = 0.1,  height = 0),
        dpi = 150
      ) + 
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                        labels = tools_levels[tools_levels %in% input$tool_overlap]) +
      scale_color_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
                         labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
      scale_x_discrete(labels = function(x) {
        x <- gsub("-", "-\n", x)
        x <- gsub(" ", "\n", x)
        x}) + 
      scale_shape_manual(values = c("no" = 16, "yes" = 17)) +
      facet_grid(tool_ref ~ new_level, scales = "free_x") +
      scale_y_continuous(limits = c(-0.1,1.1), 
                         breaks = seq(0, 1, length.out = 3),
                         labels = scales::label_number()) + 
      ylab("CSNO") + 
      xlab("") + 
      theme(
        legend.position = "none",
        legend.text = element_text(size = general_size ),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(15, 15, 15, 15, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        panel.background = element_rect(colour = "black", fill = NA),
        strip.text = element_blank(),
        strip.background = element_blank())
  })
  
  output$overlap_csno <- renderPlot({plot_overlap_csno()})
  output$overlap_csno_overview <- renderPlot({plot_overlap_csno()})
  
  plot_overlap_csno_summary  <- reactive({
    df <- recall_base()
    req(nrow(df) > 0)
    
    df %>%
      filter(!is.na(recall)) %>% 
      ungroup() %>% 
      group_by(tool_ref, new_level) %>% 
      summarise(fnr = median(fnr), .groups = "drop") %>%
      ggplot(aes(x = "Class medians", fill = tool_ref, y = fnr)) +
      geom_violin() +
      ggrastr::rasterise(
        geom_jitter(size = 1.5, width = 0.4, height = 0),
        dpi = 150
      ) + 
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                        labels = tools_levels[tools_levels %in% input$tool_overlap]) +
      facet_grid(tool_ref ~ ., scales = "free_x") +
      scale_y_continuous(limits = c(-0.2,1.2), 
                         breaks = seq(0, 1, length.out = 3),
                         labels = scales::label_number()) +
      scale_x_discrete(labels = function(x) {
        x <- gsub("-", "-\n", x)
        x <- gsub(" ", "\n", x)
        x}) + 
      ylab("") + 
      xlab("") + 
      theme(
        legend.position = "none",
        legend.text = element_text(size = general_size ),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(15, 15, 15, 15, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_blank(),
        panel.background = element_rect(colour = "black", fill = NA),
        strip.text = element_text(size = general_size, face = "bold"),
        strip.background = element_blank())
  })
  
  output$overlap_csno_summary <- renderPlot({plot_overlap_csno_summary()})
  output$overlap_csno_summary_overview <- renderPlot({plot_overlap_csno_summary()})
  
  plot_overlap_legend <- reactive({
    tools_order <- match(input$tool_overlap_calc, tools_levels)
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]
    
    df <- recall_base()
    req(nrow(df) > 0)
    
    grid.draw(g_legend(df %>%
                         ggplot(aes(x = recall + fnr, y = recall + fnr)) +
                         ggrastr::rasterise(
                           geom_jitter(aes(color = tool_comp, shape = texture), size = 2.5),
                           dpi = 150
                           ) + 
                         scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
                                           labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
                         scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
                         guides(
                           fill = guide_legend(
                             override.aes = list(shape = shape_tools, fill = pal_10_q[tools_levels %in% input$tool_overlap_calc])), shape = "none") +
                         labs(fill = "") +
                         theme(
                           legend.position = "bottom",
                           legend.text = element_text(size = general_size ),
                           panel.border = element_blank(),
                           panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank(),
                           plot.margin = margin(0, 0, 0, 0, unit = "pt"),
                           legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
                           legend.margin = margin(0, 0, 0, 0, unit = "pt"),
                           title = element_text(size = general_size + 2, face = "bold"),
                           axis.title = element_text(size = general_size + 1, face = "bold"),
                           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
                           axis.text.y = element_text(size = general_size),
                           panel.background = element_rect(colour = "black", fill = NA),
                           strip.text = element_blank(),
                           strip.background = element_blank(),
                           legend.key = element_blank(), 
                           legend.background = element_blank())))
  })
  
  output$overlap_legend <- renderPlot({plot_overlap_legend()})
  output$plot_cstc_legend <- renderPlot({plot_overlap_legend()}) # Reuses the logic
  output$plot_csno_legend <- renderPlot({plot_overlap_legend()}) # Reuses the logic
  
  
  # Pan-/Core-resistome proportion 
  pan_reactive <- reactive({
    thr_name <- get_thr_name(input$threshold_pan_core_id)
    df <- pan_prepped[[thr_name]]
    req(!is.null(df))
    df
  })
  
  sumcore_reactive <- reactive({
    thr_name <- get_thr_name(input$threshold_pan_core_id)
    sum_core_adjust(
      core_prepped[[thr_name]], 
      as.numeric(input$threshold_samples), 
      input$threshold_proportion
    )
  })
  
  output$plot_pan_core_proportion <- renderPlot({
    alluvial_pan_core_env(
      sumcore = sumcore_reactive() %>%
        filter(tool %in% input$tool_pan_core,
               habitat %in% input$single_environment_pan_core) %>%
        mutate(
          tool = factor(as.character(tool), 
                        levels = tools_levels[tools_levels %in% input$tool_pan_core])),
      pan     = pan_reactive() %>%
        filter(tool %in% input$tool_pan_core,
               habitat %in% input$single_environment_pan_core) %>%
        mutate(
          tool = factor(as.character(tool), 
                        levels = tools_levels[tools_levels %in% input$tool_pan_core])),
      h       = input$single_environment_pan_core,
      tools   = input$tool_pan_core,
      pal     = pal_10_complete,
      general_size = general_size,
      data_list$levels_unigenes) 
  })
}