#' @title Plot Coefficients from unitregTMB Models
#' 
#' @description 
#' Creates a coefficient plot (also known as a forest plot) to visualize point estimates 
#' and confidence intervals for one or more fitted \code{unitregTMB} models. 
#' It allows for easy model comparison, component selection (mean, precision, or inflation), 
#' and highly customizable \code{ggplot2} aesthetics.
#' 
#' @param ... One or more fitted model objects of class \code{unitregTMB}, or a list of such objects.
#' @param conf_level Numeric. The confidence level for the error bars. Default is \code{0.95}.
#' @param intercept Logical. Should the intercept(s) be included in the plot? Default is \code{TRUE}.
#' @param component Character vector specifying which model components to plot. 
#'        Options are \code{"all"}, \code{"mu"} (default), \code{"phi"}, \code{"p0"}, or \code{"p1"}.
#' @param ncol Integer. Number of columns for the faceted layout when \code{facet = TRUE}. Default is \code{2}.
#' @param labels Character vector. Custom labels for the y-axis (coefficient names). Default is \code{NULL}.
#' @param parse_labels Logical. If \code{TRUE}, parses \code{labels} as mathematical expressions. Default is \code{FALSE}.
#' @param math_labels Logical. If \code{TRUE}, automatically generates mathematical expressions for the coefficients 
#'        (e.g., beta^phi). Default is \code{FALSE}.
#' @param facet Logical. If \code{TRUE} (default), displays models in separate panels using \code{patchwork} and \code{ggh4x}. 
#'        If \code{FALSE}, overlays models in a single plot using position dodge.
#' @param palette Character string specifying the color palette. Options include \code{"default"}, 
#'        \code{"grey"}, \code{"viridis"}, \code{"brewer"}, or a \code{ggsci} palette (e.g., \code{"ggsci::npg"}).
#' @param point_size Numeric. Size of the point estimates. Default is \code{2}.
#' @param strip_text_size Numeric. Font size for the facet strip labels. Default is \code{11}.
#' @param axis_title_size Numeric. Font size for the axis titles. Default is \code{11}.
#' @param axis_text_x_size Numeric. Font size for the x-axis text. Default is \code{10}.
#' @param axis_text_y_size Numeric. Font size for the y-axis text. Default is \code{12}.
#' @param vline_color Character. Color for the vertical zero-reference line. Default is \code{"gray40"}.
#' @param point_color Character. Color for the points when \code{facet = TRUE}. Default is \code{"#0072B2"}.
#' @param errorbar_color Character. Color for the error bars when \code{facet = TRUE}. Default is \code{"gray40"}.
#' @param title Character. Main title for the plot. Default is \code{NULL}.
#' @param subtitle Character. Subtitle for the plot. Automatically generated based on \code{conf_level} if missing.
#' @param xlab Character. Label for the x-axis. Default is \code{"Estimate Value"}.
#' @param ylab Character. Label for the y-axis. Default is \code{"Coefficient"}.
#' 
#' @return A \code{ggplot} object containing the coefficient plot.
#' 
#' @details 
#' When \code{facet = TRUE}, the function requires the \code{patchwork}, \code{grid}, and \code{ggh4x} 
#' packages to arrange the panels seamlessly. When \code{facet = FALSE}, models are plotted together 
#' in a single panel and distinguished by color and shape.
#' 
#' @examples
#' \donttest{
#' # Assuming 'da' is your dataset:
#' # fit1 <- unitregTMB(Y ~ educ + refill, data = da, family = vasicek())
#' # fit2 <- unitregTMB(Y ~ educ + refill, data = da, family = kumaraswamy())
#' #
#' # # 1. Plot the location components (mu) for a single model
#' # plot_coef(fit1, component = "mu")
#' #
#' # # 2. Compare two models overlaying them in the same panel (facet = FALSE)
#' # plot_coef(fit1, fit2, component = "mu", facet = FALSE, palette = "Set1")
#' #
#' # # 3. Plot all components of a model with mathematical labels and a custom palette
#' # plot_coef(fit1, component = "all", math_labels = TRUE, palette = "ggsci::npg")
#' }
#' 
#' @importFrom stats qnorm
#' @export
plot_coef <- function(..., conf_level = 0.95, intercept = TRUE, component = "mu", 
                      ncol = 2, 
                      labels = NULL, 
                      parse_labels = FALSE, 
                      math_labels = FALSE, 
                      facet = TRUE, 
                      palette = "default",
                      point_size = 2,
                      strip_text_size = 11,
                      axis_title_size = 11,
                      axis_text_x_size = 10, 
                      axis_text_y_size = 12, 
                      vline_color = "gray40",
                      point_color = "#0072B2",
                      errorbar_color = "gray40",
                      title = NULL,
                      subtitle,
                      xlab = "Estimate Value",
                      ylab = "Coefficient") {
  
  if (missing(subtitle)) {
    subtitle <- paste0(conf_level * 100, "% Confidence Level")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The 'ggplot2' package is required.")
  if (facet && !requireNamespace("patchwork", quietly = TRUE)) stop("The 'patchwork' package is required for facet mode.")
  if (facet && !requireNamespace("grid", quietly = TRUE)) stop("The 'grid' package is required for facet mode.")
  if (facet && !requireNamespace("ggh4x", quietly = TRUE)) stop("The 'ggh4x' package is required for facet mode.")
  
  if (is.character(palette) && length(palette) == 1 && startsWith(palette, "ggsci::")) {
    if (!requireNamespace("ggsci", quietly = TRUE)) {
      stop("Package 'ggsci' is required to use this palette. Please install it with install.packages('ggsci').", call. = FALSE)
    }
  }
  
  models <- list(...)
  if (length(models) == 1 && is.list(models[[1]]) && !inherits(models[[1]], "unitregTMB")) {
    models <- models[[1]]
  }
  
  if (!all(sapply(models, inherits, "unitregTMB"))) {
    stop("All inputs must be 'unitregTMB' model objects.")
  }
  
  plot_data_list <- lapply(models, function(m) {
    sd_rep <- if (!is.null(m$sd_report)) m$sd_report else m$sdreport
    if (is.null(sd_rep)) return(NULL)
    summ <- summary(sd_rep, "fixed")
    
    fam_name <- if (is.list(m$family) && !is.null(m$family$name)) {
      m$family$name
    } else if (!is.null(m$family)) {
      as.character(m$family)[1]
    } else {
      "Unknown"
    }
    
    comp_prefixes <- c()
    if ("all" %in% component) {
      comp_prefixes <- c("beta_mu", "beta_phi", "beta_p0", "beta_p1")
    } else {
      if ("mu" %in% component) comp_prefixes <- c(comp_prefixes, "beta_mu")
      if ("phi" %in% component) comp_prefixes <- c(comp_prefixes, "beta_phi")
      if ("p0" %in% component) comp_prefixes <- c(comp_prefixes, "beta_p0")
      if ("p1" %in% component) comp_prefixes <- c(comp_prefixes, "beta_p1")
    }
    
    df_list <- list()
    for (pfx in comp_prefixes) {
      idx <- grep(paste0("^", pfx), rownames(summ))
      if (length(idx) > 0) {
        est <- summ[idx, 1]
        se  <- summ[idx, 2]
        mat_name <- switch(pfx, "beta_mu" = "X_mu", "beta_phi" = "X_phi", "beta_p0" = "X_p0", "beta_p1" = "X_p1")
        t_names <- rep(pfx, length(est))
        if (!is.null(m$obj$env$data[[mat_name]])) {
          r_names <- colnames(m$obj$env$data[[mat_name]])
          if (!is.null(r_names) && length(r_names) == length(est)) {
            t_names <- r_names
          } else {
            t_names <- paste0(pfx, "_", 1:length(est))
          }
        }
        if (pfx != "beta_mu") {
          pfx_short <- sub("beta_", "", pfx)
          t_names <- paste0(pfx_short, "_", t_names)
        }
        df_list[[pfx]] <- data.frame(
          model_raw_name = fam_name,
          tau = if(is.null(m$tau)) NA else m$tau,
          coef_name = trimws(t_names),
          estimate = est,
          std_error = se,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(df_list) > 0) do.call(rbind, df_list) else NULL
  })
  
  plot_data <- do.call(rbind, plot_data_list)
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    message("No coefficients to plot in the provided objects.")
    return(invisible(NULL))
  }
  
  orig_names_all <- unique(plot_data$coef_name)
  
  if (!is.null(labels)) {
    label_function <- if (parse_labels) parse(text = labels) else labels
  } else if (math_labels) {
    math_exprs <- character(length(orig_names_all))
    for(i in seq_along(orig_names_all)) {
      nm <- orig_names_all[i]
      if (grepl("^phi_", nm)) {
        idx <- sum(grepl("^phi_", orig_names_all[1:i])) - 1
        math_exprs[i] <- sprintf("beta[%d]^phi", idx)
      } else if (grepl("^p0_", nm)) {
        idx <- sum(grepl("^p0_", orig_names_all[1:i])) - 1
        math_exprs[i] <- sprintf("beta[%d]^{p0}", idx)
      } else if (grepl("^p1_", nm)) {
        idx <- sum(grepl("^p1_", orig_names_all[1:i])) - 1
        math_exprs[i] <- sprintf("beta[%d]^{p1}", idx)
      } else {
        idx <- sum(!grepl("^(phi|p0|p1)_", orig_names_all[1:i])) - 1
        math_exprs[i] <- sprintf("beta[%d]", idx)
      }
    }
    names(math_exprs) <- orig_names_all
    label_function <- function(x) parse(text = math_exprs[x])
  } else {
    label_function <- ggplot2::waiver()
  }
  
  if (!intercept) {
    plot_data <- plot_data[!grepl("(Intercept)", plot_data$coef_name, fixed = TRUE), ]
  }
  
  capitalize <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))
  plot_data$model_name <- mapply(function(name, tau_val) {
    if (grepl("\\(", name)) return(name) 
    dist_name <- capitalize(strsplit(as.character(name), "_")[[1]][1])
    
    if (!is.null(tau_val) && !is.na(tau_val) && is.numeric(tau_val)) {
      parameter <- paste0("tau = ", tau_val)
    } else {
      parameter <- strsplit(as.character(name), "_")[[1]][2]
    }
    if(is.na(parameter) || is.null(parameter)) return(dist_name)
    paste0(dist_name, " (", parameter, ")")
  }, plot_data$model_raw_name, plot_data$tau)
  
  z_value <- qnorm((1 + conf_level) / 2)
  plot_data$ci_lower <- plot_data$estimate - z_value * plot_data$std_error
  plot_data$ci_upper <- plot_data$estimate + z_value * plot_data$std_error
  
  all_coef_names <- rev(unique(plot_data$coef_name))
  plot_data$coef_name <- factor(plot_data$coef_name, levels = all_coef_names)
  model_order <- unique(plot_data$model_name)
  plot_data$model_name <- factor(plot_data$model_name, levels = model_order)
  
  if (facet) {
    plots_list_data <- split(plot_data, plot_data$model_name)
    individual_plots <- lapply(plots_list_data, function(current_data) {
      p <- ggplot2::ggplot(current_data, ggplot2::aes(y = coef_name, x = estimate)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = vline_color, linewidth = 0.8) +
        ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), color = errorbar_color, linewidth = 1.2) +
        ggplot2::geom_point(size = point_size, color = point_color) +
        ggh4x::facet_wrap2(~ model_name, strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = "grey90"))) +
        ggplot2::scale_y_discrete(limits = all_coef_names, labels = label_function, drop = FALSE) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          strip.text = ggplot2::element_text(face = "bold", size = strip_text_size),
          axis.text.x = ggplot2::element_text(size = axis_text_x_size), 
          axis.text.y = ggplot2::element_text(size = axis_text_y_size)  
        )
      return(p)
    })
    
    total_plots <- length(individual_plots)
    total_rows <- ceiling(total_plots / ncol)
    for (i in seq_along(individual_plots)) {
      is_first_col <- (i - 1) %% ncol == 0
      current_row <- (i - 1) %/% ncol + 1
      is_bottom_row <- current_row == total_rows || (i + ncol > total_plots)
      if (!is_first_col) individual_plots[[i]] <- individual_plots[[i]] + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
      if (!is_bottom_row) individual_plots[[i]] <- individual_plots[[i]] + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }
    
    main_grid <- patchwork::wrap_plots(individual_plots, ncol = ncol)
    
    y_label_grob <- grid::textGrob(ylab, rot = 90, gp = grid::gpar(fontsize = axis_title_size))
    x_label_grob <- grid::textGrob(xlab, gp = grid::gpar(fontsize = axis_title_size))
    
    top_row <- patchwork::wrap_elements(y_label_grob) + main_grid + patchwork::plot_layout(widths = c(0.05, 1))
    bottom_row <- patchwork::plot_spacer() + patchwork::wrap_elements(x_label_grob) + patchwork::plot_layout(widths = c(0.05, 1))
    
    final_layout <- top_row / bottom_row + patchwork::plot_layout(heights = c(1, 0.05))
    
    has_title <- !is.null(title) && title != ""
    has_subtitle <- !is.null(subtitle) && subtitle != ""
    
    if (has_title || has_subtitle) {
      final_plot <- final_layout +
        patchwork::plot_annotation(title = title, subtitle = subtitle) & 
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
          plot.subtitle = ggplot2::element_text(size = 12, hjust = 0)
        )
    } else {
      final_plot <- final_layout
    }
    
  } else {
    pd <- ggplot2::position_dodge(width = 0.6)
    
    final_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(y = coef_name, x = estimate, color = model_name, shape = model_name)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = vline_color, linewidth = 0.8) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), position = pd, linewidth = 1) +
      ggplot2::geom_point(size = point_size, position = pd, stroke = 1.2) +
      ggplot2::scale_y_discrete(limits = all_coef_names, labels = label_function) +
      ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = ylab, color = "Model", shape = "Model") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16), 
        legend.position = "top",
        axis.title = ggplot2::element_text(size = axis_title_size),
        axis.text.x = ggplot2::element_text(size = axis_text_x_size), 
        axis.text.y = ggplot2::element_text(size = axis_text_y_size)  
      )
    
    num_models <- length(levels(plot_data$model_name))
    
    if (is.character(palette) && length(palette) > 1) {
      final_plot <- final_plot + ggplot2::scale_color_manual(values = palette)
    } else if (is.character(palette) && length(palette) == 1 && startsWith(palette, "ggsci::")) {
      ggsci_palette_name <- sub("ggsci::", "", palette)
      scale_function_name <- paste0("scale_color_", ggsci_palette_name)
      
      if (exists(scale_function_name, where = asNamespace("ggsci"), mode = "function")) {
        ggsci_scale_function <- get(scale_function_name, asNamespace("ggsci"))
        final_plot <- final_plot + ggsci_scale_function()
      } else {
        warning("ggsci palette '", ggsci_palette_name, "' not found. Using default colors.", call. = FALSE)
      }
      
    } else {
      predefined_palettes <- c("default", "grey", "viridis", "brewer")
      if (palette %in% predefined_palettes) {
        final_plot <- final_plot + switch(
          palette,
          "grey"    = list(ggplot2::scale_color_grey(start = 0.0, end = 0.6), 
                           ggplot2::scale_shape_manual(values = 1:num_models)),
          "viridis" = ggplot2::scale_color_viridis_d(),
          "brewer"  = ggplot2::scale_color_brewer(palette = "Set2"),
          NULL 
        )
      } else {
        warning("Palette '", palette, "' not recognized. Using default ggplot2 colors.\n", call. = FALSE)
      }
    }
  }
  
  return(final_plot)
}
