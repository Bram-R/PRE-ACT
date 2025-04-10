f_stm_diagram <- function(state_labels, transitions) {
  #' Create a State Diagram
  #'
  #' This function generates a state diagram using DiagrammeR, with customizable state labels
  #' and transitions. It ensures accurate arrow placements, including self-loops.
  #'
  #' @param state_labels Named list. A list of state labels with keys as state names
  #'   and values as positions (e.g., list("state1" = c("Event free", "2,2!"))).
  #' @param transitions Data frame. A data frame defining the transitions between states.
  #'   Columns should include `from`, `to`, and optionally `label` for edge attributes.
  #' @return An HTML widget displaying the state diagram.
  #' @examples
  #' labels <- list(
  #'   state1 = c("Event free", "2,2!"),
  #'   state2 = c("Locoregional recurrence", "0,0!"),
  #'   state3 = c("Distant metastasis", "4,0!"),
  #'   state4 = c("Death", "2,-2!")
  #' )
  #' transitions <- data.frame(
  #' from = c("state1", "state1", "state1", "state1",
  #'          "state2", "state2", "state2",
  #'          "state3", "state3",
  #'          "state4"),
  #' to = c("state1", "state2", "state3", "state4",
  #'        "state2", "state3", "state4",
  #'        "state3", "state4",
  #'        "state4"),
  #' label = rep("", 10) 
  #' )
  #' f_stm_diagram(labels, transitions)
  #' @importFrom DiagrammeR grViz
  #' @export
  
  source("f_stm_diagram_helper.R") # Load helper functions
  
  validate_state_labels(state_labels)
  validate_transitions(transitions, state_labels)
  
  nodes <- generate_nodes(state_labels)
  edges <- generate_edges(transitions, state_labels)
  layout_constraints <- "graph [layout = neato];"
  
  diagram_code <- paste(
    "digraph state_diagram {",
    nodes,
    edges,
    layout_constraints,
    "}",
    sep = "\n"
  )
  
  DiagrammeR::grViz(diagram_code)
}