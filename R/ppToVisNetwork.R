#' Converts a PanelPRO compatible data frame into a visNetwork data frame
#'
#' @param ped a pedigree that contains \code{ID, Sex, MotherID, FatherID, isProband, CurAge}
#' It is assumed that in this data.frame, Sex = 1 if male, Sex = 0 if female and Sex = NA if unkonwn
#' Affliction is coded in "isAffXX" columns, where XX are short names of cancers
#' The age of diagnosis is correspondingly coded in "AgeXX" columns
#' @param title a string for the title on the plot
#' @import dplyr
#' @import tidyr
#' @importFrom rlang :=
#' @import visNetwork
#' @export
toVisNetwork <- function(ped, title = "Your Pedigree") {
  # Get the counselee(s)
  counseleeIDs <- ped %>% filter(isProband == 1) %>% select(ID) %>%
    unlist() %>% as.vector()

  # Get the cancer affections
  cancersInPedigree <- substring(colnames(ped[, grepl("isAff", names(ped)), drop = FALSE]), 6)

  # Get cancer columns
  pedCancers <- ped %>% select(ID, starts_with("isAff") | starts_with("Age"))

  for (cancer in cancersInPedigree) {
    pedCancers <- pedCancers %>% mutate("c.{cancer}" := if_else(get(paste0("isAff", cancer)) == 1,
                                                                cancer,
                                                                NA_character_)) %>%
      mutate("ca.{cancer}" := if_else(get(paste0("isAff", cancer)) == 1 & is.na(get(paste0("Age", cancer))),
                                      paste0(cancer, ": age missing"),
                                      # paste0(cancer, ": ", get(paste0("ca.", cancer)))))
                                      NA_character_)) %>%
      mutate("ca.{cancer}" := if_else(get(paste0("isAff", cancer)) == 1 & !is.na(get(paste0("Age", cancer))),
                                      paste0(cancer, ": ", as.character(get(paste0("Age", cancer)))),
                                      get(paste0("ca.", cancer))))

  }

  # Concatenate into strings
  pedCancers <- pedCancers %>% unite(group, starts_with("c."), remove = TRUE,
                                     sep = ", ", na.rm = TRUE) %>%
    unite(ages, starts_with("ca."), remove = TRUE, sep = "<br>", na.rm = TRUE) %>%
    mutate(group = if_else(group == "", "unaffected", group)) %>%
    select(ID, group, ages)

  # visNetwork requires an id column and from/to columns
  nodes <- data.frame(id = ped$ID,
                      shape = ifelse(ped$Sex, "box", "ellipse"),
                      label = paste0(ifelse(ped$isDead, "D: ", "C: "),
                                     ped$CurAge),
                      stringsAsFactors = FALSE) %>%
    left_join(pedCancers, by = c("id" = "ID")) %>%
    mutate(group = if_else(id %in% counseleeIDs, paste0(group, ", counselee"),
                           group)) %>%
    mutate(title = paste0("ID: ", id, "<br>")) %>%
    mutate(title = paste0(title, ifelse(ped$isDead, "Dead at: ", "Current age: "),
                          ped$CurAge, "<br>", ages))

  connections <- ped %>% select(ID, MotherID, FatherID) %>%
    filter(!is.na(MotherID) & !is.na(FatherID))

  edges <- bind_rows(
    data.frame(from = connections$MotherID, to = connections$ID),
    data.frame(from = connections$FatherID, to = connections$ID)
  )

  visNetwork(nodes, edges, main = title) %>%
    visEdges(arrows = "to",
             color = list(color = "gray")) %>%
    visOptions(selectedBy = list(variable = "group", multiple = TRUE),
               nodesIdSelection = TRUE,
               highlightNearest = list(enabled = TRUE,
                                       degree = 2,
                                       hover = TRUE),
               manipulation = list(enabled = FALSE)) %>%
    visLegend(addNodes = list(list(label = "D: Death age", shape = "icon",
                                   icon = list(code = "f007", color = "gray",
                                               size = 20)),
                              list(label = "C: Current age", shape = "icon",
                                   icon = list(code = "f007", color = "gray",
                                               size = 20))),
              useGroups = TRUE,
              ncol = 2) %>%
    visHierarchicalLayout(sortMethod = "directed")
}

# Global variables
utils::globalVariables(c("isProband", "ID", "MotherID", "FatherID",
                         "group", "ages"))
