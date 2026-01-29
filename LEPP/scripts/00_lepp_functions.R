###############################################################################
# LEPP: Low Expression as Predictors of Phenotype
# Function library
#
# This script defines all statistical models, bootstrapping utilities,
# and result summarization helpers used in the LEPP DKO vs SKO analysis.
#
# No data loading, plotting, or side effects should occur in this file.
###############################################################################

## ============================================================================
## Utility helpers
## ============================================================================

# Extract estimate and standard error for a specific model term
# ----------------------------------------------------------------
# tidy_df   : broom.mixed::tidy() output
# term_name : name of coefficient to extract
#
# Returns:
#   list(estimate = numeric, se = numeric) or NULL if missing
extract_term_est_se <- function(tidy_df, term_name) {
    
    # Ensure required columns exist
    if (!all(c("term", "estimate", "std.error") %in% colnames(tidy_df))) {
        return(NULL)
    }
    
    # Subset to requested term
    row <- tidy_df[tidy_df$term == term_name, , drop = FALSE]
    
    # Return NULL if term is absent
    if (nrow(row) == 0) return(NULL)
    
    # Coerce to numeric and return
    list(
        estimate = as.numeric(row$estimate[1]),
        se       = as.numeric(row$std.error[1])
    )
}

## ============================================================================
## Formula construction
## ============================================================================

# Build a regression formula that automatically drops invalid covariates
# -----------------------------------------------------------------------
# df            : data.frame containing model variables
# response      : name of response variable (string)
# fixed_effects : vector of candidate fixed-effect covariates
# random_effect : optional random effect (string)
#
# Returns:
#   stats::formula object safe for glmmTMB
build_safe_formula <- function(df, response, fixed_effects, random_effect = NULL) {
    
    # Helper: determine if a covariate has usable variance
    covariate_is_valid <- function(x) {
        
        # Factor / character: must have >1 observed level
        if (is.factor(x) || is.character(x)) {
            length(unique(x[!is.na(x)])) > 1
            
            # Numeric: must have non-zero variance
        } else if (is.numeric(x)) {
            stats::var(x, na.rm = TRUE) > 0
            
            # Otherwise invalid
        } else {
            FALSE
        }
    }
    
    # Retain only fixed effects that pass validity checks
    valid_fixed <- fixed_effects[
        vapply(
            fixed_effects,
            function(v) covariate_is_valid(df[[v]]),
            logical(1)
        )
    ]
    
    # Determine whether random effect is usable
    use_random <- FALSE
    if (!is.null(random_effect) && random_effect %in% colnames(df)) {
        use_random <- length(unique(df[[random_effect]])) > 1
    }
    
    # Assemble RHS terms
    rhs_terms <- valid_fixed
    
    # Append random intercept if valid
    if (use_random) {
        rhs_terms <- c(rhs_terms, paste0("(1 | ", random_effect, ")"))
    }
    
    # Construct RHS string
    rhs <- if (length(rhs_terms) == 0) "1" else paste(rhs_terms, collapse = " + ")
    
    # Return formula
    as.formula(paste(response, "~", rhs))
}

## ============================================================================
## Single-gene model
## ============================================================================

# Fit a single-gene low-expression model
# -------------------------------------
# df    : analysis data.frame
# geneA : gene symbol (string)
#
# Returns:
#   list(geneA, geneB = NA, tidy, model) or NULL on failure
fit_gene_single_soft_low <- function(df, geneA) {
    
    # Depth-normalized log-rate
    df$A_lr <- log1p(df[[geneA]]) - log(df$nCount_RNA)
    
    # Define lowness: higher values = lower expression
    df$A_low <- -df$A_lr
    
    # Z-score lowness
    df$single_low_z <- scale(df$A_low)[, 1]
    
    # Build safe model formula
    f <- build_safe_formula(
        df = df,
        response = "phenotype_score",
        fixed_effects = c("single_low_z", "subset", "tissue"),
        random_effect = "patient_id"
    )
    
    # Fit model, fail gracefully
    fit <- tryCatch(
        glmmTMB::glmmTMB(formula = f, family = gaussian(), data = df),
        error = function(e) NULL
    )
    
    # Exit if model failed
    if (is.null(fit)) return(NULL)
    
    # Return structured output
    list(
        geneA = geneA,
        geneB = NA_character_,
        tidy  = broom.mixed::tidy(fit),
        model = fit
    )
}

## ============================================================================
## Pairwise (DKO-style) model
## ============================================================================

# Fit a paired low-expression interaction model
# ---------------------------------------------
# df    : analysis data.frame
# geneA : first gene
# geneB : second gene
#
# Returns:
#   list(geneA, geneB, tidy, model) or NULL on failure
fit_gene_pair_soft_both_low <- function(df, geneA, geneB) {
    
    # Depth-normalized log-rates
    df$A_lr <- log1p(df[[geneA]]) - log(df$nCount_RNA)
    df$B_lr <- log1p(df[[geneB]]) - log(df$nCount_RNA)
    
    # Convert to lowness scale
    df$A_low <- -df$A_lr
    df$B_low <- -df$B_lr
    
    # BOTH-low penalty: dominated by the lower of the two
    df$pair_low_penalty <- pmin(df$A_low, df$B_low)
    
    # Flip sign so positive beta = stronger phenotype at low expression
    df$pair_low_z <- -scale(df$pair_low_penalty)[, 1]
    
    # Include main effects for interpretability
    df$A_low_z <- scale(df$A_low)[, 1]
    df$B_low_z <- scale(df$B_low)[, 1]
    
    # Build safe model formula
    f <- build_safe_formula(
        df = df,
        response = "phenotype_score",
        fixed_effects = c("A_low_z", "B_low_z", "pair_low_z", "subset", "tissue"),
        random_effect = "patient_id"
    )
    
    # Fit model
    fit <- tryCatch(
        glmmTMB::glmmTMB(formula = f, family = gaussian(), data = df),
        error = function(e) NULL
    )
    
    # Exit on failure
    if (is.null(fit)) return(NULL)
    
    # Return structured output
    list(
        geneA = geneA,
        geneB = geneB,
        tidy  = broom.mixed::tidy(fit),
        model = fit
    )
}

## ============================================================================
## Result summarization
## ============================================================================

# Summarize pairwise model outputs into a tidy table
summarize_pair_results <- function(res) {
    
    # Drop failed models
    res <- Filter(Negate(is.null), res)
    if (length(res) == 0) return(tibble::tibble())
    
    dplyr::bind_rows(
        lapply(res, function(x) {
            x$tidy |>
                dplyr::mutate(
                    geneA = x$geneA,
                    geneB = x$geneB
                )
        })
    ) |>
        dplyr::filter(term %in% c("A_low_z", "B_low_z", "pair_low_z")) |>
        dplyr::mutate(
            effect = dplyr::case_when(
                term == "A_low_z"    ~ "GeneA_low",
                term == "B_low_z"    ~ "GeneB_low",
                term == "pair_low_z" ~ "BothLow"
            )
        )
}

# Summarize single-gene model outputs
summarize_single_results <- function(res) {
    
    # Drop failed models
    res <- Filter(Negate(is.null), res)
    if (length(res) == 0) return(tibble::tibble())
    
    dplyr::bind_rows(
        lapply(res, function(x) {
            x$tidy |>
                dplyr::mutate(
                    geneA     = x$geneA,
                    geneB     = x$geneB,
                    component = "magnitude"
                )
        })
    )
}

## ============================================================================
## Parametric bootstrap: pair vs best single
## ============================================================================

# Fast parametric bootstrap comparing pair effect vs best single
# --------------------------------------------------------------
# mu_single  : vector of single-gene estimates
# se_single  : vector of single-gene SEs
# pair_tbl   : data.frame with mu_pair and se_pair
#
# Returns:
#   pair_tbl augmented with delta statistics
bootstrap_pair_vs_best_single_parametric <- function(
        mu_single, se_single, pair_tbl,
        n_boot = 1000, seed = 1, chunk_size = 250
) {
    
    set.seed(seed)
    
    n_s <- length(mu_single)
    n_p <- nrow(pair_tbl)
    
    # Early exit if no usable inputs
    if (n_s == 0 || n_p == 0) {
        return(
            pair_tbl |>
                dplyr::mutate(
                    delta_median = NA_real_,
                    delta_ci_low = NA_real_,
                    delta_ci_high = NA_real_,
                    prob_pair_gt_single = NA_real_,
                    delta_mean = NA_real_,
                    delta_sd   = NA_real_
                )
        )
    }
    
    # Split bootstrap into chunks for parallelization
    starts <- seq(1, n_boot, by = chunk_size)
    chunks <- lapply(starts, function(s) s:min(s + chunk_size - 1, n_boot))
    
    # Run bootstrap chunks in parallel
    chunk_delta_list <- future.apply::future_lapply(
        chunks,
        function(idx) {
            
            B <- length(idx)
            
            # Sample single-gene effects
            single_draws <- matrix(rnorm(B * n_s), B, n_s)
            single_draws <- sweep(single_draws, 2, se_single, `*`)
            single_draws <- sweep(single_draws, 2, mu_single, `+`)
            best_single  <- apply(single_draws, 1, max, na.rm = TRUE)
            
            # Sample pair effects
            pair_draws <- matrix(rnorm(B * n_p), B, n_p)
            pair_draws <- sweep(pair_draws, 2, pair_tbl$se_pair, `*`)
            pair_draws <- sweep(pair_draws, 2, pair_tbl$mu_pair, `+`)
            
            # Delta: pair – best single
            pair_draws - best_single
        },
        future.seed = TRUE
    )
    
    # Combine all bootstrap samples
    delta_mat <- do.call(rbind, chunk_delta_list)
    
    # Summarize bootstrap distribution
    pair_tbl |>
        dplyr::mutate(
            delta_median = apply(delta_mat, 2, median, na.rm = TRUE),
            delta_ci_low = apply(delta_mat, 2, quantile, 0.025, na.rm = TRUE),
            delta_ci_high = apply(delta_mat, 2, quantile, 0.975, na.rm = TRUE),
            prob_pair_gt_single = apply(delta_mat, 2, function(x) mean(x > 0, na.rm = TRUE)),
            delta_mean = apply(delta_mat, 2, mean, na.rm = TRUE),
            delta_sd   = apply(delta_mat, 2, sd,   na.rm = TRUE)
        )
}

###############################################################################
# End of function library
###############################################################################
