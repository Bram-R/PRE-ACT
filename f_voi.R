f_make_model_voi <- function(model_fn = f_model, par_fn = f_input, setting = n_setting) {
  #' Create a model function compatible with `voi::evppi_mc()`
  #'
  #' Creates a wrapper around a state-transition model that accepts a single
  #' parameter data frame (e.g., `f_model(params)`) and converts it to the
  #' interface required by `voi::evppi_mc()`, where each uncertain parameter is
  #' supplied as a separate function argument.
  #'
  #' @param model_fn Function evaluating a single parameter set (e.g. `f_model()`).
  #' @param par_fn Function generating model inputs (e.g. `f_input()`).
  #' @param setting Integer specifying the model setting passed to `par_fn()`.
  #'
  #' @return A function compatible with `voi::evppi_mc()`.
  #'
  #' @export
  
  # Obtain parameter names
  v_par_names <- names(
    par_fn(
      n_sim = 1,
      setting = setting
    )
  )
  
  # Generate wrapper function
  eval(
    parse(
      text = sprintf(
        "
        function(%s) {

          params <- data.frame(
            %s,
            check.names = FALSE
          )

          res <- model_fn(params)

          rbind(
            Effect = res[(n_treatments + 1):(2 * n_treatments)],
            Cost   = res[1:n_treatments]
          )
        }
        ",
        paste(v_par_names, collapse = ", "),
        paste(
          sprintf("%s = %s", v_par_names, v_par_names),
          collapse = ",\n            "
        )
      )
    )
  )
}