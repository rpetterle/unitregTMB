#' Plot Coefficients from unitregTMB Models
#'
#' @description
#' Creates coefficient plots with Wald confidence intervals for one or more
#' fitted \code{unitregTMB} models. Fixed effects from the \code{mu},
#' \code{phi}, \code{p0}, and \code{p1} distributional predictors can be
#' displayed separately or together.
#'
#' @param ... One or more fitted objects of class \code{unitregTMB}, or a
#'   single list containing such objects.
#' @param conf_level Numeric scalar in \eqn{(0,1)} giving the confidence level
#'   for the Wald intervals.
#' @param intercept Logical; include intercept coefficients.
#' @param component Character vector selecting \code{"mu"}, \code{"phi"},
#'   \code{"p0"}, or \code{"p1"}. Use \code{"all"} to include all active
#'   components.
#' @param ncol Positive integer giving the number of columns when
#'   \code{facet = TRUE}.
#' @param labels Optional character vector used to relabel coefficients. A
#'   named vector may be used to map coefficient names to labels.
#' @param parse_labels Logical; parse custom labels as plotmath expressions.
#' @param math_labels Logical; automatically use mathematical labels
#'   \eqn{\beta}, \eqn{\psi}, \eqn{\gamma}, and \eqn{\delta} for coefficients
#'   in the \code{mu}, \code{phi}, \code{p0}, and \code{p1} predictors,
#'   respectively.
#' @param facet Logical; if \code{TRUE}, models are shown in separate panels.
#'   If \code{FALSE}, models are overlaid in one panel.
#' @param palette Color specification used when \code{facet = FALSE}. It may be
#'   \code{"default"}, \code{"grey"}, \code{"viridis"}, \code{"brewer"}, a
#'   character vector of colors, or a \code{ggsci} scale specified as
#'   \code{"ggsci::npg"}, for example.
#' @param point_size Numeric point size.
#' @param strip_text_size Numeric font size for facet strip labels.
#' @param axis_title_size Numeric font size for axis titles.
#' @param axis_text_x_size Numeric font size for x-axis labels.
#' @param axis_text_y_size Numeric font size for y-axis labels.
#' @param vline_color Color of the zero-reference line.
#' @param point_color Point color used when \code{facet = TRUE}.
#' @param errorbar_color Confidence-interval color used when
#'   \code{facet = TRUE}.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle. If omitted, the confidence level is
#'   displayed.
#' @param xlab Label for the x axis.
#' @param ylab Label for the y axis.
#'
#' @details
#' Confidence intervals are Wald intervals computed from the fixed-effect
#' estimates and standard errors stored in the fitted \code{unitregTMB}
#' object. Therefore, the model must have been fitted with standard-error
#' calculation enabled.
#'
#' When several models have the same family label, unique suffixes are added
#' automatically so that their coefficients are not combined accidentally.
#'
#' @return A \code{ggplot} object when \code{facet = FALSE}. When
#'   \code{facet = TRUE}, a \code{patchwork} object containing the faceted
#'   coefficient plots is returned.
#'
#' @examples
#' \donttest{
#' data("wqi", package = "unitregTMB")
#'
#' fit1 <- unitregTMB(
#'   y ~ location,
#'   data = wqi,
#'   family = beta_fam(model_for = "mean")
#' )
#'
#' fit2 <- unitregTMB(
#'   y ~ location,
#'   data = wqi,
#'   family = unitgamma(model_for = "mean")
#' )
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_coef(fit1, component = "mu", facet = FALSE)
#'   plot_coef(fit1, fit2, component = "mu", facet = FALSE)
#' }
#' }
#'
#' @author Ricardo Rasmussen Petterle
#' @importFrom stats qnorm
#' @export
plot_coef <- function(..., conf_level = 0.95, intercept = TRUE,
                      component = "mu", ncol = 2, labels = NULL,
                      parse_labels = FALSE, math_labels = FALSE,
                      facet = TRUE, palette = "default", point_size = 2,
                      strip_text_size = 11, axis_title_size = 11,
                      axis_text_x_size = 10, axis_text_y_size = 12,
                      vline_color = "gray40", point_color = "#0072B2",
                      errorbar_color = "gray40", title = NULL, subtitle,
                      xlab = "Estimate Value", ylab = "Coefficient") {

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).", call. = FALSE)
  }

  if (length(intercept) != 1L || !is.logical(intercept) || is.na(intercept)) {
    stop("'intercept' must be TRUE or FALSE.", call. = FALSE)
  }

  valid_components <- c("all", "mu", "phi", "p0", "p1")
  if (!is.character(component) || !length(component) ||
      any(!component %in% valid_components)) {
    stop(
      "'component' must contain only 'all', 'mu', 'phi', 'p0', or 'p1'.",
      call. = FALSE
    )
  }

  if ("all" %in% component) {
    component <- c("mu", "phi", "p0", "p1")
  } else {
    component <- unique(component)
  }

  if (!is.numeric(ncol) || length(ncol) != 1L || !is.finite(ncol) ||
      ncol < 1 || ncol != as.integer(ncol)) {
    stop("'ncol' must be a positive integer.", call. = FALSE)
  }
  ncol <- as.integer(ncol)

  for (x in c("parse_labels", "math_labels", "facet")) {
    value <- get(x, inherits = FALSE)
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop("'", x, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }

  if (missing(subtitle)) {
    subtitle <- paste0(format(100 * conf_level, trim = TRUE), "% Confidence Level")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required.", call. = FALSE)
  }

  if (facet && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("The 'patchwork' package is required when 'facet = TRUE'.", call. = FALSE)
  }

  if (facet && !requireNamespace("ggh4x", quietly = TRUE)) {
    stop("The 'ggh4x' package is required when 'facet = TRUE'.", call. = FALSE)
  }

  if (is.character(palette) && length(palette) == 1L &&
      startsWith(palette, "ggsci::") &&
      !requireNamespace("ggsci", quietly = TRUE)) {
    stop("The 'ggsci' package is required for this palette.", call. = FALSE)
  }

  models <- list(...)
  if (length(models) == 1L && is.list(models[[1L]]) &&
      !inherits(models[[1L]], "unitregTMB")) {
    models <- models[[1L]]
  }

  if (!length(models)) {
    stop("At least one fitted 'unitregTMB' model must be supplied.", call. = FALSE)
  }

  ok <- vapply(models, inherits, logical(1), what = "unitregTMB")
  if (!all(ok)) {
    stop("All inputs must inherit from 'unitregTMB'.", call. = FALSE)
  }

  model_base_label <- function(m) {
    fam <- if (is.list(m$family) && !is.null(m$family$name)) {
      as.character(m$family$name)[1L]
    } else if (!is.null(m$family) && length(m$family)) {
      as.character(m$family)[1L]
    } else {
      "unitregTMB"
    }

    if (!is.null(m$tau) && length(m$tau) == 1L && is.finite(m$tau) &&
        !grepl("tau\\s*=", fam)) {
      fam <- paste0(fam, " (tau = ", format(m$tau, trim = TRUE), ")")
    }

    fam
  }

  model_labels <- make.unique(
    vapply(models, model_base_label, character(1)),
    sep = " #"
  )

  plot_data_list <- lapply(seq_along(models), function(i) {
    m <- models[[i]]
    rows <- vector("list", length(component))
    names(rows) <- component

    for (comp in component) {
      active <- switch(
        comp,
        mu = TRUE,
        phi = TRUE,
        p0 = isTRUE(m$has_p0),
        p1 = isTRUE(m$has_p1)
      )

      if (!active) next

      est <- m$model_coef[[comp]]
      if (is.null(est) || !length(est)) next

      se <- m$model_std[[paste0("std.", comp)]]
      if (is.null(se)) se <- rep(NA_real_, length(est))

      if (length(se) != length(est)) {
        stop(
          "Coefficient and standard-error lengths differ for component '",
          comp, "' in model '", model_labels[i], "'.",
          call. = FALSE
        )
      }

      raw_names <- names(est)
      if (is.null(raw_names) || length(raw_names) != length(est)) {
        raw_names <- paste0("coef", seq_along(est))
      }

      display_names <- if (comp == "mu") {
        raw_names
      } else {
        paste0(comp, "_", raw_names)
      }

      rows[[comp]] <- data.frame(
        model_name = model_labels[i],
        component = comp,
        coef_index = seq_along(est) - 1L,
        raw_name = raw_names,
        coef_name = display_names,
        estimate = as.numeric(est),
        std_error = as.numeric(se),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (!length(rows)) return(NULL)
    do.call(rbind, rows)
  })

  plot_data_list <- plot_data_list[!vapply(plot_data_list, is.null, logical(1))]
  if (!length(plot_data_list)) {
    stop("No requested active coefficients are available to plot.", call. = FALSE)
  }

  plot_data <- do.call(rbind, plot_data_list)
  rownames(plot_data) <- NULL

  if (!intercept) {
    plot_data <- plot_data[plot_data$raw_name != "(Intercept)", , drop = FALSE]
  }

  if (!nrow(plot_data)) {
    stop("No coefficients remain after applying 'intercept = FALSE'.", call. = FALSE)
  }

  if (any(!is.finite(plot_data$estimate)) ||
      any(!is.finite(plot_data$std_error))) {
    stop(
      "Finite fixed-effect estimates and standard errors are required. ",
      "Refit with standard-error calculation enabled.",
      call. = FALSE
    )
  }

  z_value <- stats::qnorm((1 + conf_level) / 2)
  plot_data$ci_lower <- plot_data$estimate - z_value * plot_data$std_error
  plot_data$ci_upper <- plot_data$estimate + z_value * plot_data$std_error

  coef_levels <- rev(unique(plot_data$coef_name))
  plot_data$coef_name <- factor(plot_data$coef_name, levels = coef_levels)
  plot_data$model_name <- factor(plot_data$model_name, levels = model_labels)

  if (!is.null(labels)) {
    if (!is.character(labels)) {
      stop("'labels' must be a character vector.", call. = FALSE)
    }

    if (is.null(names(labels))) {
      target <- unique(as.character(plot_data$coef_name))
      target <- target[!is.na(target)]
      if (length(labels) != length(target)) {
        stop(
          "Unnamed 'labels' must have one entry for each plotted coefficient.",
          call. = FALSE
        )
      }
      names(labels) <- target
    }

    label_map <- setNames(coef_levels, coef_levels)
    hit <- intersect(names(labels), names(label_map))
    label_map[hit] <- labels[hit]

    if (parse_labels) {
      label_function <- function(x) parse(text = unname(label_map[x]))
    } else {
      label_function <- function(x) unname(label_map[x])
    }
  } else if (math_labels) {
    meta <- plot_data[
      !duplicated(as.character(plot_data$coef_name)),
      c("coef_name", "component", "coef_index"),
      drop = FALSE
    ]

    symbols <- c(mu = "beta", phi = "psi", p0 = "gamma", p1 = "delta")
    math_map <- setNames(
      paste0(symbols[meta$component], "[", meta$coef_index, "]"),
      as.character(meta$coef_name)
    )

    label_function <- function(x) parse(text = unname(math_map[x]))
  } else {
    label_function <- ggplot2::waiver()
  }

  if (facet) {
    plots_list_data <- split(plot_data, plot_data$model_name, drop = TRUE)

    individual_plots <- lapply(plots_list_data, function(current_data) {
      ggplot2::ggplot(
        current_data,
        ggplot2::aes(y = coef_name, x = estimate)
      ) +
        ggplot2::geom_vline(
          xintercept = 0, linetype = "dashed",
          color = vline_color, linewidth = 0.8
        ) +
        ggplot2::geom_linerange(
          ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
          color = errorbar_color, linewidth = 1.2
        ) +
        ggplot2::geom_point(size = point_size, color = point_color) +
        ggh4x::facet_wrap2(
          ~ model_name,
          strip = ggh4x::strip_themed(
            background_x = ggh4x::elem_list_rect(fill = "grey90")
          )
        ) +
        ggplot2::scale_y_discrete(
          limits = coef_levels,
          labels = label_function,
          drop = FALSE
        ) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          strip.text = ggplot2::element_text(
            face = "bold", size = strip_text_size
          ),
          axis.text.x = ggplot2::element_text(size = axis_text_x_size),
          axis.text.y = ggplot2::element_text(size = axis_text_y_size)
        )
    })

    total_plots <- length(individual_plots)
    total_rows <- ceiling(total_plots / ncol)

    for (i in seq_along(individual_plots)) {
      is_first_col <- (i - 1L) %% ncol == 0L
      current_row <- (i - 1L) %/% ncol + 1L
      is_bottom_row <- current_row == total_rows || (i + ncol > total_plots)

      if (!is_first_col) {
        individual_plots[[i]] <- individual_plots[[i]] +
          ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )
      }

      if (!is_bottom_row) {
        individual_plots[[i]] <- individual_plots[[i]] +
          ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank()
          )
      }
    }

    main_grid <- patchwork::wrap_plots(individual_plots, ncol = ncol)
    y_label_grob <- grid::textGrob(
      ylab, rot = 90, gp = grid::gpar(fontsize = axis_title_size)
    )
    x_label_grob <- grid::textGrob(
      xlab, gp = grid::gpar(fontsize = axis_title_size)
    )

    top_row <- patchwork::wrap_elements(y_label_grob) +
      main_grid +
      patchwork::plot_layout(widths = c(0.05, 1))

    bottom_row <- patchwork::plot_spacer() +
      patchwork::wrap_elements(x_label_grob) +
      patchwork::plot_layout(widths = c(0.05, 1))

    final_plot <- top_row / bottom_row +
      patchwork::plot_layout(heights = c(1, 0.05))

    has_title <- !is.null(title) && length(title) == 1L && nzchar(title)
    has_subtitle <- !is.null(subtitle) && length(subtitle) == 1L && nzchar(subtitle)

    if (has_title || has_subtitle) {
      final_plot <- final_plot +
        patchwork::plot_annotation(title = title, subtitle = subtitle) &
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold", size = 16, hjust = 0
          ),
          plot.subtitle = ggplot2::element_text(size = 12, hjust = 0)
        )
    }
  } else {
    pd <- ggplot2::position_dodge(width = 0.6)

    final_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        y = coef_name, x = estimate,
        color = model_name, shape = model_name
      )
    ) +
      ggplot2::geom_vline(
        xintercept = 0, linetype = "dashed",
        color = vline_color, linewidth = 0.8
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
        position = pd, linewidth = 1
      ) +
      ggplot2::geom_point(
        size = point_size, position = pd, stroke = 1.2
      ) +
      ggplot2::scale_y_discrete(
        limits = coef_levels, labels = label_function
      ) +
      ggplot2::labs(
        title = title, subtitle = subtitle,
        x = xlab, y = ylab, color = "Model", shape = "Model"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 16),
        legend.position = "top",
        axis.title = ggplot2::element_text(size = axis_title_size),
        axis.text.x = ggplot2::element_text(size = axis_text_x_size),
        axis.text.y = ggplot2::element_text(size = axis_text_y_size)
      )

    num_models <- length(unique(plot_data$model_name))

    if (is.character(palette) && length(palette) > 1L) {
      final_plot <- final_plot +
        ggplot2::scale_color_manual(values = palette)
    } else if (is.character(palette) && length(palette) == 1L &&
               startsWith(palette, "ggsci::")) {
      palette_name <- sub("^ggsci::", "", palette)
      scale_name <- paste0("scale_color_", palette_name)

      if (exists(scale_name, where = asNamespace("ggsci"), mode = "function")) {
        scale_fun <- get(scale_name, envir = asNamespace("ggsci"))
        final_plot <- final_plot + scale_fun()
      } else {
        warning(
          "ggsci palette '", palette_name,
          "' was not found; using default colors.",
          call. = FALSE
        )
      }
    } else if (is.character(palette) && length(palette) == 1L) {
      if (identical(palette, "grey")) {
        final_plot <- final_plot +
          ggplot2::scale_color_grey(start = 0, end = 0.6) +
          ggplot2::scale_shape_manual(values = seq_len(num_models))
      } else if (identical(palette, "viridis")) {
        final_plot <- final_plot + ggplot2::scale_color_viridis_d()
      } else if (identical(palette, "brewer")) {
        final_plot <- final_plot +
          ggplot2::scale_color_brewer(palette = "Set2")
      } else if (!identical(palette, "default")) {
        warning(
          "Palette '", palette,
          "' was not recognized; using default ggplot2 colors.",
          call. = FALSE
        )
      }
    }
  }

  final_plot
}
