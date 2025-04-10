#### f_stm_diagram helper functions ----

# check function using: docstring(validate_state_labels)
# check function using: docstring(validate_transitions)
# check function using: docstring(generate_nodes)
# check function using: docstring(generate_edges)
# check function using: docstring(get_self_loop_style)

validate_state_labels <- function(state_labels) {
  #' Validate State Labels
  #'
  #' Ensures the state labels input is correctly formatted.
  #'
  #' @param state_labels Named list of state names and their labels/positions.
  #' @return NULL. Throws an error if validation fails.
  if (!is.list(state_labels) || length(state_labels) < 2) {
    stop("state_labels must be a named list with at least two states.")
  }
  
  for (label in state_labels) {
    if (length(label) != 2 || !is.character(label)) {
      stop("Each state label must be a character vector of length 2: c(label, position).")
    }
  }
}

validate_transitions <- function(transitions, state_labels) {
  #' Validate Transitions
  #'
  #' Ensures the transitions data frame is correctly formatted.
  #'
  #' @param transitions Data frame. Transition definitions.
  #' @param state_labels Named list. State labels for validation.
  #' @return NULL. Throws an error if validation fails.
  required_columns <- c("from", "to")
  if (!all(required_columns %in% colnames(transitions))) {
    stop("Transitions data frame must have at least 'from' and 'to' columns.")
  }
  
  if (!all(transitions$from %in% names(state_labels)) ||
      !all(transitions$to %in% names(state_labels))) {
    stop("All transitions must refer to valid states in state_labels.")
  }
}

generate_nodes <- function(state_labels) {
  #' Generate Node Definitions
  #'
  #' Creates the node definitions for the state diagram.
  #'
  #' @param state_labels Named list of state names and their labels/positions.
  #' @return A string of node definitions in Graphviz syntax.
  paste(
    sapply(names(state_labels), function(state) {
      sprintf(
        "%s [label = '%s', pos = '%s'];",
        state, state_labels[[state]][1], state_labels[[state]][2]
      )
    }),
    collapse = "\n"
  )
}

generate_edges <- function(transitions, state_labels) {
  #' Generate Edge Definitions
  #'
  #' Creates the edge definitions for the state diagram based on transitions.
  #'
  #' @param transitions Data frame. Defines the transitions between states.
  #' @param state_labels Named list of state names.
  #' @return A string of edge definitions in Graphviz syntax.
  paste(
    apply(transitions, 1, function(row) {
      from <- row["from"]
      to <- row["to"]
      label <- ifelse(!is.na(row["label"]), row["label"], "")
      edge_style <- if (from == to) {
        get_self_loop_style(from)
      } else {
        "fontsize = 10, arrowsize = 0.75, penwidth = 0.75"
      }
      
      sprintf(
        "%s -> %s [label = '%s', %s];",
        from, to, label, edge_style
      )
    }),
    collapse = "\n"
  )
}

get_self_loop_style <- function(state) {
  #' Get Self-Loop Style
  #'
  #' Provides the style for self-looping transitions based on node position.
  #'
  #' @param state Character. State name.
  #' @return A string specifying Graphviz edge attributes for self-loops.
  ports <- list(
    state1 = "n", state2 = "w", state3 = "e", state4 = "s"
  )
  port <- ports[[state]]
  sprintf("fontsize = 10, arrowsize = 0.75, penwidth = 0.75, headport = %s, tailport = %s", port, port)
}
