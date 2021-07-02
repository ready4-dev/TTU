#' Predict from shareable
#' @description predict_from_shareable_mdl() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict from shareable model. The function returns New data (a double vector).
#' @param model_mdl Model (a model)
#' @param data_tb Data (a tibble)
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: 'response'
#' @param sd_dbl Standard deviation (a double vector)
#' @param deterministic_1L_lgl Deterministic (a logical vector of length one), Default: T
#' @return New data (a double vector)
#' @rdname predict_from_shareable_mdl
#' @export 
#' @importFrom rlang exec
predict_from_shareable_mdl <- function (model_mdl, data_tb, predn_type_1L_chr = "response", 
    sd_dbl, deterministic_1L_lgl = T) 
{
    pred_fn <- ifelse(inherits(model_mdl, "betareg"), predict_shrble_betareg, 
        ifelse(inherits(model_mdl, "glm"), predict_shrble_glm, 
            predict_shrble_lm))
    args_ls <- list(object = model_mdl, newdata = data_tb, type = predn_type_1L_chr, 
        sd_1L_dbl = max(0, rnorm(1, mean = sd_dbl[1], sd = ifelse(deterministic_1L_lgl, 
            0, sd_dbl[2]))))
    new_data_dbl <- rlang::exec(pred_fn, !!!args_ls)
    return(new_data_dbl)
}
#' Predict shrble betareg
#' @description predict_shrble_betareg() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict shrble betareg. The function is called for its side effects and does not return a value.
#' @param object PARAM_DESCRIPTION
#' @param newdata PARAM_DESCRIPTION, Default: NULL
#' @param type PARAM_DESCRIPTION, Default: c("response", "link", "precision", "variance", "quantile")
#' @param na.action PARAM_DESCRIPTION, Default: na.pass
#' @param at PARAM_DESCRIPTION, Default: 0.5
#' @param sd_1L_dbl Standard deviation (a double vector of length one)
#' @param ... Additional arguments
#' @return NULL
#' @rdname predict_shrble_betareg
#' @export 

predict_shrble_betareg <- function (object, newdata = NULL, type = c("response", "link", 
    "precision", "variance", "quantile"), na.action = na.pass, 
    at = 0.5, sd_1L_dbl, ...) 
{
    type <- match.arg(type)
    if (type == "quantile") {
        qfun <- function(at, mu, phi) {
            rval <- sapply(at, function(p) qbeta(p, mu * phi, 
                (1 - mu) * phi))
            if (length(at) > 1L) {
                if (NCOL(rval) == 1L) 
                  rval <- matrix(rval, ncol = length(at), dimnames = list(unique(names(rval)), 
                    NULL))
                colnames(rval) <- paste("q_", at, sep = "")
            }
            else {
                rval <- drop(rval)
            }
            rval
        }
    }
    if (missing(newdata)) {
        rval <- switch(type, response = {
            object$fitted.values
        }, link = {
            object$link$mean$linkfun(object$fitted.values)
        }, precision = {
            gamma <- object$coefficients$precision
            z <- if (is.null(object$x)) model.matrix(object, 
                model = "precision") else object$x$precision
            offset <- if (is.null(object$offset$precision)) rep.int(0, 
                NROW(z)) else object$offset$precision
            object$link$precision$linkinv(drop(z %*% gamma + 
                offset))
        }, variance = {
            gamma <- object$coefficients$precision
            z <- if (is.null(object$x)) model.matrix(object, 
                model = "precision") else object$x$precision
            offset <- if (is.null(object$offset$precision)) rep.int(0, 
                NROW(z)) else object$offset$precision
            phi <- object$link$precision$linkinv(drop(z %*% gamma + 
                offset))
            object$fitted.values * (1 - object$fitted.values)/(1 + 
                phi)
        }, quantile = {
            mu <- object$fitted.values
            gamma <- object$coefficients$precision
            z <- if (is.null(object$x)) model.matrix(object, 
                model = "precision") else object$x$precision
            offset <- if (is.null(object$offset$precision)) rep.int(0, 
                NROW(z)) else object$offset$precision
            phi <- object$link$precision$linkinv(drop(z %*% gamma + 
                offset))
            qfun(at, mu, phi)
        })
        return(rval)
    }
    else {
        tnam <- switch(type, response = "mean", link = "mean", 
            precision = "precision", variance = "full", quantile = "full")
        mf <- model.frame(delete.response(object$terms[[tnam]]), 
            newdata, na.action = na.action, xlev = object$levels[[tnam]])
        newdata <- newdata[rownames(mf), , drop = FALSE]
        offset <- list(mean = rep.int(0, nrow(mf)), precision = rep.int(0, 
            nrow(mf)))
        if (type %in% c("response", "link", "variance", "quantile")) {
            X <- model.matrix(delete.response(object$terms$mean), 
                mf, contrasts = object$contrasts$mean)
            if (!is.null(object$call$offset)) 
                offset[[1L]] <- offset[[1L]] + eval(object$call$offset, 
                  newdata)
            if (!is.null(off.num <- attr(object$terms$mean, "offset"))) {
                for (j in off.num) offset[[1L]] <- offset[[1L]] + 
                  eval(attr(object$terms$mean, "variables")[[j + 
                    1L]], newdata)
            }
            X[, 1] <- X[, 1] + rnorm(length(X[, 1]), mean = 0, 
                sd = sd_1L_dbl)
        }
        if (type %in% c("precision", "variance", "quantile")) {
            Z <- model.matrix(object$terms$precision, mf, contrasts = object$contrasts$precision)
            if (!is.null(off.num <- attr(object$terms$precision, 
                "offset"))) {
                for (j in off.num) offset[[2L]] <- offset[[2L]] + 
                  eval(attr(object$terms$precision, "variables")[[j + 
                    1L]], newdata)
            }
        }
        rval <- switch(type, response = {
            object$link$mean$linkinv(drop(X %*% object$coefficients$mean + 
                offset[[1L]]))
        }, link = {
            drop(X %*% object$coefficients$mean + offset[[1L]])
        }, precision = {
            object$link$precision$linkinv(drop(Z %*% object$coefficients$precision + 
                offset[[2L]]))
        }, variance = {
            mu <- object$link$mean$linkinv(drop(X %*% object$coefficients$mean + 
                offset[[1L]]))
            phi <- object$link$precision$linkinv(drop(Z %*% object$coefficients$precision + 
                offset[[2L]]))
            mu * (1 - mu)/(1 + phi)
        }, quantile = {
            mu <- object$link$mean$linkinv(drop(X %*% object$coefficients$mean + 
                offset[[1L]]))
            phi <- object$link$precision$linkinv(drop(Z %*% object$coefficients$precision + 
                offset[[2L]]))
            qfun(at, mu, phi)
        })
        return(rval)
    }
}
#' Predict shrble glm
#' @description predict_shrble_glm() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict shrble glm. The function is called for its side effects and does not return a value.
#' @param object PARAM_DESCRIPTION
#' @param newdata PARAM_DESCRIPTION, Default: NULL
#' @param type PARAM_DESCRIPTION, Default: c("link", "response", "terms")
#' @param se.fit PARAM_DESCRIPTION, Default: FALSE
#' @param dispersion PARAM_DESCRIPTION, Default: NULL
#' @param terms PARAM_DESCRIPTION, Default: NULL
#' @param na.action PARAM_DESCRIPTION, Default: na.pass
#' @param sd_1L_dbl Standard deviation (a double vector of length one)
#' @param ... Additional arguments
#' @return NULL
#' @rdname predict_shrble_glm
#' @export 

predict_shrble_glm <- function (object, newdata = NULL, type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL, 
    na.action = na.pass, sd_1L_dbl, ...) 
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors, 
                response = object$fitted.values, terms = predict.lm(object, 
                  se.fit = se.fit, scale = 1, type = "terms", 
                  terms = terms))
            if (!is.null(na.act)) 
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predict_shrble_lm(object, newdata, se.fit, 
                scale = 1, type = if (type == "link") 
                  "response"
                else type, terms = terms, sd_1L_dbl = sd_1L_dbl, 
                na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    }
    else {
        if (inherits(object, "survreg")) 
            dispersion <- 1
        if (is.null(dispersion) || dispersion == 0) 
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predict_shrble_lm(object, newdata, se.fit, scale = residual.scale, 
            type = if (type == "link") 
                "response"
            else type, terms = terms, sd_1L_dbl = sd_1L_dbl, 
            na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}
#' Predict shrble linear model
#' @description predict_shrble_lm() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict shrble linear model. The function is called for its side effects and does not return a value.
#' @param object PARAM_DESCRIPTION
#' @param newdata PARAM_DESCRIPTION
#' @param se.fit PARAM_DESCRIPTION, Default: FALSE
#' @param scale PARAM_DESCRIPTION, Default: NULL
#' @param df Data.frame (a data.frame), Default: Inf
#' @param interval PARAM_DESCRIPTION, Default: c("none", "confidence", "prediction")
#' @param level PARAM_DESCRIPTION, Default: 0.95
#' @param type PARAM_DESCRIPTION, Default: c("response", "terms")
#' @param terms PARAM_DESCRIPTION, Default: NULL
#' @param na.action PARAM_DESCRIPTION, Default: na.pass
#' @param pred.var PARAM_DESCRIPTION, Default: res.var/weights
#' @param weights PARAM_DESCRIPTION, Default: 1
#' @param sd_1L_dbl Standard deviation (a double vector of length one)
#' @param ... Additional arguments
#' @return NULL
#' @rdname predict_shrble_lm
#' @export 

predict_shrble_lm <- function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
    interval = c("none", "confidence", "prediction"), level = 0.95, 
    type = c("response", "terms"), terms = NULL, na.action = na.pass, 
    pred.var = res.var/weights, weights = 1, sd_1L_dbl, ...) 
{
    tt <- terms(object)
    if (!inherits(object, "lm")) 
        warning("calling predict.lm(<fake-lm-object>) ...")
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action, 
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(tt, 
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset)) 
            offset <- offset + eval(object$call$offset, newdata)
        mmDone <- FALSE
        X[, 1] <- X[, 1] + rnorm(length(X[, 1]), mean = 0, sd = sd_1L_dbl)
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq_len(p)
    qr_fn <- function(x, ...) {
        if (is.null(r <- x$qr)) 
            stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
        r
    }
    piv <- if (p) 
        qr_fn(object)$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    if (!is.null(offset)) 
        predictor <- predictor + offset
    interval <- match.arg(interval)
    if (interval == "prediction") {
        if (missing(newdata)) 
            warning("predictions on current data refer to _future_ responses\n")
        if (missing(newdata) && missing(weights)) {
            w <- weights.default(object)
            if (!is.null(w)) {
                weights <- w
                warning("assuming prediction variance inversely proportional to weights used for fitting\n")
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) && 
            missing(pred.var)) 
            warning("Assuming constant prediction variance even though model fit is weighted\n")
        if (inherits(weights, "formula")) {
            if (length(weights) != 2L) 
                stop("'weights' as formula should be one-sided")
            d <- if (missing(newdata) || is.null(newdata)) 
                model.frame(object)
            else newdata
            weights <- eval(weights[[2L]], d, environment(weights))
        }
    }
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        w <- object$weights
        res.var <- if (is.null(scale)) {
            r <- object$residuals
            rss <- sum(if (is.null(w)) r^2 else r^2 * w)
            df <- object$df.residual
            rss/df
        }
        else scale^2
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w)) 
                  qr.Q(qr.lm(object))[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(qr.lm(object))[p1, 
                  p1])
                ip <- drop(XRinv^2 %*% rep(res.var, p))
            }
            else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrix(object)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept) 
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            avx <- colMeans(mm)
            termsconst <- sum(avx[piv] * beta[piv])
        }
        nterms <- length(asgn)
        if (nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))
            if (se.fit || interval != "none") {
                ip <- matrix(ncol = nterms, nrow = NROW(X))
                dimnames(ip) <- list(rownames(X), names(asgn))
                Rinv <- qr.solve(qr.R(qr.lm(object))[p1, p1])
            }
            if (hasintercept) 
                X <- sweep(X, 2L, avx, check.margin = FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0L] <- 0L
                predictor[, i] <- if (any(iipiv > 0L)) 
                  X[, iipiv, drop = FALSE] %*% beta[iipiv]
                else 0
                if (se.fit || interval != "none") 
                  ip[, i] <- if (any(iipiv > 0L)) 
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii, 
                      , drop = FALSE])^2 %*% rep.int(res.var, 
                      p)
                  else 0
            }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit) 
                  ip <- ip[, terms, drop = FALSE]
            }
        }
        else {
            predictor <- ip <- matrix(0, n, 0L)
        }
        attr(predictor, "constant") <- if (hasintercept) 
            termsconst
        else 0
    }
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip), 
            prediction = sqrt(ip + pred.var))
        if (type != "terms") {
            predictor <- cbind(predictor, predictor + hwid %o% 
                c(1, -1))
            colnames(predictor) <- c("fit", "lwr", "upr")
        }
        else {
            if (!is.null(terms)) 
                hwid <- hwid[, terms, drop = FALSE]
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    if (se.fit || interval != "none") {
        se <- sqrt(ip)
        if (type == "terms" && !is.null(terms) && !se.fit) 
            se <- se[, terms, drop = FALSE]
    }
    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
        predictor <- napredict(na.act, predictor)
        if (se.fit) 
            se <- napredict(na.act, se)
    }
    if (type == "terms" && interval != "none") {
        if (missing(newdata) && !is.null(na.act)) {
            lwr <- napredict(na.act, lwr)
            upr <- napredict(na.act, upr)
        }
        list(fit = predictor, se.fit = se, lwr = lwr, upr = upr, 
            df = df, residual.scale = sqrt(res.var))
    }
    else if (se.fit) 
        list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
    else predictor
}
#' Predict uncnstrd utility
#' @description predict_uncnstrd_utl() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict uncnstrd utility. The function returns New data (a double vector).
#' @param data_tb Data (a tibble)
#' @param model_mdl Model (a model)
#' @param new_data_is_1L_chr New data is (a character vector of length one), Default: 'Predicted'
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @param tfmn_for_bnml_1L_lgl Transformation for binomial (a logical vector of length one), Default: F
#' @param deterministic_1L_lgl Deterministic (a logical vector of length one), Default: T
#' @param family_1L_chr Family (a character vector of length one), Default: 'NA'
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param is_brms_mdl_1L_lgl Is bayesian regression models model (a logical vector of length one), Default: F
#' @param force_new_data_1L_lgl Force new data (a logical vector of length one), Default: F
#' @param sd_dbl Standard deviation (a double vector), Default: NA
#' @return New data (a double vector)
#' @rdname predict_uncnstrd_utl
#' @export 
#' @importFrom stats predict simulate rnorm sigma
#' @importFrom brms posterior_predict
#' @importFrom rlang exec
#' @importFrom enrichwith get_simulate_function enrich
predict_uncnstrd_utl <- function (data_tb, model_mdl, new_data_is_1L_chr = "Predicted", 
    predn_type_1L_chr = NULL, tfmn_for_bnml_1L_lgl = F, deterministic_1L_lgl = T, 
    family_1L_chr = NA_character_, tfmn_1L_chr = "NTF", is_brms_mdl_1L_lgl = F, 
    force_new_data_1L_lgl = F, sd_dbl = NA_real_) 
{
    if (new_data_is_1L_chr == "Predicted") 
        if (!is.na(sd_dbl[1])) {
            new_data_dbl <- predict_from_shareable_mdl(model_mdl, 
                data_tb = data_tb, predn_type_1L_chr = predn_type_1L_chr, 
                sd_dbl = sd_dbl, deterministic_1L_lgl = deterministic_1L_lgl)
        }
        else {
            new_data_dbl <- stats::predict(model_mdl, newdata = data_tb, 
                type = predn_type_1L_chr)
        }
    if (new_data_is_1L_chr == "Simulated") {
        if (is_brms_mdl_1L_lgl) {
            new_data_dbl <- brms::posterior_predict(model_mdl, 
                newdata = data_tb, nsamples = 1) %>% as.vector()
        }
        else {
            if ("betareg" %in% class(model_mdl) & !force_new_data_1L_lgl) {
                new_data_dbl <- rlang::exec(enrichwith::get_simulate_function(model_mdl), 
                  coef(enrichwith::enrich(model_mdl, with = "auxiliary functions")))
            }
            else {
                if (!tfmn_for_bnml_1L_lgl & !force_new_data_1L_lgl) {
                  new_data_dbl <- stats::simulate(model_mdl)$sim_1
                }
                else {
                  if (!is.na(sd_dbl[1])) {
                    new_data_dbl <- predict_from_shareable_mdl(model_mdl, 
                      data_tb = data_tb, predn_type_1L_chr = predn_type_1L_chr, 
                      sd_dbl = sd_dbl, deterministic_1L_lgl = deterministic_1L_lgl)
                  }
                  else {
                    new_data_dbl <- stats::predict(model_mdl, 
                      newdata = data_tb, type = predn_type_1L_chr)
                  }
                  if (!"betareg" %in% class(model_mdl)) 
                    new_data_dbl <- new_data_dbl + stats::rnorm(nrow(data_tb), 
                      0, stats::sigma(model_mdl))
                }
            }
        }
    }
    if (is.matrix(new_data_dbl)) 
        new_data_dbl <- new_data_dbl[, 1]
    new_data_dbl <- new_data_dbl %>% calculate_dpnt_var_tfmn(tfmn_1L_chr = ifelse(tfmn_for_bnml_1L_lgl & 
        new_data_is_1L_chr == "Simulated", ifelse(family_1L_chr == 
        "quasibinomial(log)", "LOG", ifelse(family_1L_chr == 
        "quasibinomial(logit)", "LOGIT", ifelse(family_1L_chr == 
        "quasibinomial(cloglog)", "CLL", "NTF"))), tfmn_1L_chr), 
        tfmn_is_outp_1L_lgl = T)
    return(new_data_dbl)
}
#' Predict utility
#' @description predict_utility() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict utility. The function returns Predicted utility (a double vector).
#' @param data_tb Data (a tibble)
#' @param tfmn_1L_chr Transformation (a character vector of length one), Default: 'NTF'
#' @param model_mdl Model (a model)
#' @param force_min_max_1L_lgl Force minimum maximum (a logical vector of length one), Default: T
#' @param force_new_data_1L_lgl Force new data (a logical vector of length one), Default: F
#' @param utl_min_val_1L_dbl Utility minimum value (a double vector of length one), Default: 0.03
#' @param impute_1L_lgl Impute (a logical vector of length one), Default: T
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @param new_data_is_1L_chr New data is (a character vector of length one), Default: 'Predicted'
#' @param predn_type_1L_chr Prediction type (a character vector of length one), Default: NULL
#' @param sd_dbl Standard deviation (a double vector), Default: NA
#' @param tfmn_for_bnml_1L_lgl Transformation for binomial (a logical vector of length one), Default: F
#' @param family_1L_chr Family (a character vector of length one), Default: 'NA'
#' @param is_brms_mdl_1L_lgl Is bayesian regression models model (a logical vector of length one), Default: T
#' @return Predicted utility (a double vector)
#' @rdname predict_utility
#' @export 
#' @importFrom rlang exec
predict_utility <- function (data_tb, tfmn_1L_chr = "NTF", model_mdl, force_min_max_1L_lgl = T, 
    force_new_data_1L_lgl = F, utl_min_val_1L_dbl = 0.03, impute_1L_lgl = T, 
    utl_cls_fn = NULL, new_data_is_1L_chr = "Predicted", predn_type_1L_chr = NULL, 
    sd_dbl = NA_real_, tfmn_for_bnml_1L_lgl = F, family_1L_chr = NA_character_, 
    is_brms_mdl_1L_lgl = T) 
{
    predd_utl_dbl <- predict_uncnstrd_utl(data_tb = data_tb, 
        model_mdl = model_mdl, new_data_is_1L_chr = new_data_is_1L_chr, 
        predn_type_1L_chr = predn_type_1L_chr, sd_dbl = sd_dbl, 
        tfmn_for_bnml_1L_lgl = tfmn_for_bnml_1L_lgl, family_1L_chr = family_1L_chr, 
        tfmn_1L_chr = tfmn_1L_chr, force_new_data_1L_lgl = force_new_data_1L_lgl, 
        is_brms_mdl_1L_lgl = is_brms_mdl_1L_lgl)
    if (impute_1L_lgl) 
        predd_utl_dbl[which(is.na(predd_utl_dbl))] <- predd_utl_dbl %>% 
            na.omit() %>% mean()
    if (force_min_max_1L_lgl) {
        predd_utl_dbl[which(predd_utl_dbl > 1)] <- 1
        predd_utl_dbl[which(predd_utl_dbl < utl_min_val_1L_dbl)] <- utl_min_val_1L_dbl
    }
    if (!is.null(utl_cls_fn)) {
        predd_utl_dbl <- predd_utl_dbl %>% rlang::exec(.fn = utl_cls_fn)
    }
    return(predd_utl_dbl)
}
#' Predict utility from k10
#' @description predict_utl_from_k10() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict utility from k10. The function is called for its side effects and does not return a value.
#' @param k10_1L_dbl K10 (a double vector of length one)
#' @param b0_aqol_mdl_1L_dbl B0 Assessment of Quality of Life model (a double vector of length one), Default: 0.204665
#' @param b1_aqol_mdl_1L_dbl B1 Assessment of Quality of Life model (a double vector of length one), Default: -3.617134
#' @param b0_eq5d_mdl_1L_dbl B0 eq5d model (a double vector of length one), Default: 0.8644649
#' @param b1_eq5d_mdl_1L_dbl B1 eq5d model (a double vector of length one), Default: -2.926161
#' @param aqol_error_1L_dbl Assessment of Quality of Life error (a double vector of length one), Default: 0
#' @param eq5d_error_1L_dbl Eq5d error (a double vector of length one), Default: 0
#' @return NA ()
#' @rdname predict_utl_from_k10
#' @export 

predict_utl_from_k10 <- function (k10_1L_dbl, b0_aqol_mdl_1L_dbl = 0.204665, b1_aqol_mdl_1L_dbl = -3.617134, 
    b0_eq5d_mdl_1L_dbl = 0.8644649, b1_eq5d_mdl_1L_dbl = -2.926161, 
    aqol_error_1L_dbl = 0, eq5d_error_1L_dbl = 0) 
{
    meanaqol8dutility <- exp(b0_aqol_mdl_1L_dbl + b1_aqol_mdl_1L_dbl * 
        k10_1L_dbl * 0.01) + aqol_error_1L_dbl
    if (is.na(meanaqol8dutility)) 
        stop("Mean utility calculation is returning NAs")
    meaneq5dutility <- b0_eq5d_mdl_1L_dbl + b1_eq5d_mdl_1L_dbl * 
        (k10_1L_dbl * 0.01)^2 + eq5d_error_1L_dbl
    if (is.na(meaneq5dutility)) 
        stop("Mean EQ5D utility calculation is returning NAs")
    return(c(meanaqol8dutility, meaneq5dutility))
}
