gibbs_sldax = function(m = 100, burn = 0, docs, w, K = 2L,
                       model = c("lda", "slda", "sldax", "slda_logit", "sldax_logit"),
                       y = NULL, x = NULL, interaction_xcol = -1L,
                       alpha_ = 0.1, gamma_ = 1.01,
                       mu0 = NULL, sigma0 = NULL, a0 = NULL, b0 = NULL,
                       eta_start = NULL, constrain_eta = TRUE,
                       proposal_sd = NULL,
                       verbose = FALSE, display_progress = FALSE) {

  # Check model and convert to integer codes
  if (is.character(model)) {
    if (length(model) == 1) {
      model = model
    } else {
      model = model[1]
    }
    if (!model %in% c("lda", "slda", "sldax", "slda_logit", "sldax_logit")) {
      print(model)
      stop("'model' not recognized")
    }
  } else {
    print(model)
    stop("'model' not recognized")
  }
  model = switch(model,
                 lda = 1L,
                 slda = 2L,
                 sldax = 3L,
                 slda_logit = 4L,
                 sldax_logit = 5L)

  if (model == 1) {
    if (!is.null(y) || !is.null(x) || !is.null(mu0) || !is.null(sigma0) ||
        !is.null(a0) || !is.null(b0) || !is.null(eta_start) ||
        !is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'lda' model")
    }
  }
  if (model == 2) {
    if (!is.null(x) || !is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'slda' model")
    }
  }
  if (model == 3) {
    if (!is.null(proposal_sd)) {
      stop("Invalid parameter supplied for 'sldax' model")
    }
  }
  if (model == 4) {
    if (!is.null(x) || !is.null(a0) || !is.null(b0)) {
      stop("Invalid parameter supplied for 'slda_logit' model")
    }
  }
  if (model == 5) {
    if (!is.null(a0) || !is.null(b0)) {
      stop("Invalid parameter supplied for 'sldax_logit' model")
    }
  }

  # Default priors for supervised models
  if (model %in% c(2, 4)) {
    if (is.null(mu0)) mu0 = rep(0, K)
    if (is.null(sigma0)) {
      if (model == 2) sigma0 = diag(1, K)
      if (model == 4) sigma0 = diag(2.5 ^ 2, K) # Roughly equivalent to recommended Cauchy prior scales of 2.5 from Gelman paper
    }
    if (is.null(eta_start)) eta_start = seq(1, -1, length.out = K)
  }
  if (model %in% c(3, 5)) {
    if (interaction_xcol <= 0) {
      p = ncol(x)
    } else {
      p = ncol(x) + K - 1
    }
    if (is.null(mu0)) mu0 = rep(0, p + K)
    if (is.null(sigma0)) {
      if (model == 3) sigma0 = diag(1, p + K)
      if (model == 5) sigma0 = diag(2.5 ^ 2, p + K) # Roughly equivalent to recommended Cauchy prior scales of 2.5 from Gelman paper
    }
    if (is.null(eta_start)) eta_start = seq(1, -1, length.out = p + K)
  }
  if (model %in% c(2, 3)) {
    if (is.null(a0)) a0 = 0.001
    if (is.null(b0)) b0 = 0.001
  }
  if (model == 4) {
    if (is.null(proposal_sd)) proposal_sd = rep(2.38, K)
  }
  if (model == 5) {
    if (interaction_xcol <= 0) {
      p = ncol(x)
    } else {
      p = ncol(x) + K - 1
    }
    if (is.null(proposal_sd)) proposal_sd = rep(2.38, p + K)
  }

  # Set up unused defaults for C++ call
  if (model == 1) {
    y = -1L
    x = matrix(0, 1, 1)
    mu0 = 0
    sigma0 = matrix(0, 1, 1)
    a0 = -1
    b0 = -1
    eta_start = 0
    proposal_sd = 0
  }
  if (model == 2) {
    x = matrix(0, 1, 1)
    proposal_sd = 0
  }
  if (model == 3) {
    proposal_sd = 0
  }
  if (model == 4) {
    x = matrix(0, 1, 1)
    a0 = -1
    b0 = -1
  }
  if (model == 5) {
    a0 = -1
    b0 = -1
  }

  # Check m
  if (is.numeric(m)) {
    if (!is.na(m) & !is.infinite(m)) {
      if (m %% 1 == 0) {
        m = as.integer(m)
      }
    } else {
      print(m)
      stop("'m' is not an integer")
    }
  } else {
    print(m)
    stop("'m' is not an integer")
  }

  # Check burn
  if (is.numeric(burn)) {
    if (!is.na(burn) & !is.infinite(burn)) {
      if (burn %% 1 == 0) {
        burn = as.integer(burn)
      }
    } else {
      print(burn)
      stop("'burn' is not an integer")
    }
  } else {
    print(burn)
    stop("'burn' is not an integer")
  }

  # Check number of topics
  if (is.numeric(K)) {
    if (!is.na(K) & !is.infinite(K)) {
      if (K %% 1 == 0) {
        K = as.integer(K)
      }
    } else {
      print(K)
      stop("'K' is not an integer")
    }
  } else {
    print(K)
    stop("'K' is not an integer")
  }

  gibbs_sldax_cpp(docs, w, m, burn, K, model, y, x, mu0, sigma0, a0, b0,
                  eta_start, proposal_sd, interaction_xcol, alpha_, gamma_,
                  constrain_eta, verbose, display_progress)
}
