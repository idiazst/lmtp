# Supplied splits use wide-row indices; cluster IDs must stay together.
validate_folds <- function(folds, n, id) {
  if (!length(folds)) stop("`folds` must contain at least one fold.")
  for (fold in folds) {
    for (name in c("training_set", "validation_set")) {
      rows <- fold[[name]]
      if (!is.numeric(rows) || !length(rows) || anyNA(rows) ||
          any(rows != as.integer(rows)) || any(rows < 1 | rows > n) || anyDuplicated(rows)) {
        stop("Each fold needs nonempty, unique `training_set` and `validation_set` row indices.")
      }
    }
    if (length(folds) > 1) {
      if (!setequal(fold$training_set, setdiff(seq_len(n), fold$validation_set))) {
        stop("Each training set must be the complement of its validation set.")
      }
      if (any(id[fold$training_set] %in% id[fold$validation_set])) {
        stop("Cluster IDs cannot occur in both training and validation sets.")
      }
    } else if (!setequal(fold$training_set, seq_len(n)) ||
               !setequal(fold$validation_set, seq_len(n))) {
      stop("A single fold must use all rows for both training and validation.")
    }
  }
  rows <- unlist(lapply(folds, `[[`, "validation_set"), use.names = FALSE)
  if (anyDuplicated(rows) || !setequal(rows, seq_len(n))) {
    stop("Validation sets must partition all rows exactly once.")
  }
  folds
}

validate_learner_folds <- function(learner_folds, folds, id) {
  if (is.null(learner_folds)) return(NULL)
  if (!is.list(learner_folds) || length(learner_folds) != length(folds)) {
    stop("`learner_folds` must contain one fold-label vector per outer training set.")
  }
  for (j in seq_along(folds)) {
    labels <- learner_folds[[j]]
    train_id <- id[folds[[j]]$training_set]
    if (!is.numeric(labels) || length(labels) != length(train_id) || anyNA(labels) ||
        any(labels < 1 | labels != as.integer(labels))) {
      stop("Inner fold labels must be positive integers in outer training-row order.")
    }
    if (any(vapply(split(labels, train_id), function(x) length(unique(x)) != 1L, logical(1)))) {
      stop("All rows from a cluster must have the same inner fold label.")
    }
  }
  learner_folds
}

learner_split <- function(task, fold, ids, default) {
  if (is.null(task$learner_folds)) return(default)
  train_id <- task$id[task$folds[[fold]]$training_set]
  labels <- task$learner_folds[[fold]][match(ids, train_id)]
  if (anyNA(labels)) stop("Learner CV contains IDs outside the outer training set.")
  list(fold_id = labels)
}

prepare_curve_splits <- function(task, control) {
  if (!is.null(task$learner_folds)) return(invisible(task))
  V <- max(c(control$.learners_trt_folds, control$.learners_outcome_folds, 2L))
  task$learner_folds <- lapply(task$folds, function(fold) {
    rows <- fold$training_set
    if (V > length(unique(task$id[rows]))) {
      stop("The number of learner CV folds exceeds the number of training clusters; reduce the learner fold counts.")
    }
    splits <- make_folds(task$natural[rows, , drop = FALSE], V = V, cluster_ids = task$id[rows])
    labels <- integer(length(rows))
    for (j in seq_along(splits)) labels[splits[[j]]$validation_set] <- j
    labels
  })
  invisible(task)
}

# Training pseudo-outcomes must not use nuisance fits trained on the outer
# validation fold. Refit their weights with inner cross-fitting inside T_j.
curve_training_weights <- function(task, fold, ratios, sporadic_weights, learners, mtp, control) {
  rows <- task$folds[[fold]]$training_set
  if (length(task$folds) == 1L) {
    return(list(ratios = ratios[rows, , drop = FALSE],
                sporadic = sporadic_weights[rows, , drop = FALSE]))
  }
  inner <- task$clone(deep = TRUE)
  inner$natural <- task$natural[rows, , drop = FALSE]
  inner$shifted <- task$shifted[rows, , drop = FALSE]
  inner$id <- task$id[rows]
  inner$n <- length(rows)
  inner$weights <- task$weights[rows]
  labels <- task$learner_folds[[fold]]
  if (length(unique(labels)) < 2L) stop("Cross-fitted SDR needs at least two inner folds.")
  inner$folds <- lapply(unique(labels), function(label) {
    list(training_set = which(labels != label), validation_set = which(labels == label))
  })
  inner$learner_folds <- NULL
  list(ratios = cf_r(inner, learners, mtp, control, function(...) NULL)$ratios,
       sporadic = cf_sporadic(inner, learners, control, function(...) NULL)$weights)
}
