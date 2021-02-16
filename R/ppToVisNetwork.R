#' Converts a PanelPRO compatible data frame into a visNetwork data frame
#'
#' @param ped a pedigree that contains \code{ID, Sex, MotherID, FatherID, isProband, CurAge}
#' It is assumed that in this data.frame, Sex = 1 if male, Sex = 0 if female and Sex = NA if unkonwn
#' Affliction is coded in "isAffXX" columns, where XX are short names of cancers
#' The age of diagnosis is correspondingly coded in "AgeXX" columns
#' @param title a string for the title on the plot
#' @export
toVisNetwork <- function(ped, title = "Your Pedigree") {
  # Get the counselee(s)
  # Colour the counselee(s) black
  counseleeIDs <- ped %>% filter(isProband == 1) %>% select(ID)

  # visNetwork requires an id column and from/to columns
  nodes <- data.frame(id = ped$ID,
                      shape = "icon",
                      icon.code = ifelse(ped$Sex,
                                         "f183",
                                         "f182"),
                      icon.color = ifelse(ped$isProband,
                                          "black",
                                          "gray"),
                      # label = paste0("ID: ", ped$ID),
                      label = paste0(ifelse(ped$isDead, "D: ", "C: "),
                                     ped$CurAge),
                      stringsAsFactors = FALSE
                      ) %>%
    mutate(title = paste0(ifelse(ped$isDead, "Dead at: ", "Current age: "),
                          ped$CurAge))

  connections <- ped %>% select(ID, MotherID, FatherID) %>%
    filter(!is.na(MotherID) & !is.na(FatherID))

  edges <- bind_rows(
    data.frame(from = connections$MotherID, to = connections$ID),
    data.frame(from = connections$FatherID, to = connections$ID)
  )

  # Get the cancer affections



  visNetwork(nodes, edges, main = title) %>%
    visEdges(arrows = "to",
             color = list(color = "gray")) %>%
    visOptions(nodesIdSelection = FALSE,
               highlightNearest = list(enabled = TRUE,
                                       degree = 2,
                                       hover = TRUE),
               manipulation = list(enabled = FALSE)) %>%
                                   # editEdgeCols = c("label"),
                                   # editNodeCols = list(
                                   #   "text" = c("id", "label", "title"),
                                   #   "number" = c("size")),
                                   # addNodeCols = c("label", "group"))) %>%
    visLegend(addNodes = list(list(label = "counselee", shape = "icon",
                                   icon = list(code = "f007", color = "black",
                                               size = 20)),
                              list(label = "D: Death age", shape = "icon",
                                   icon = list(code = "f007", color = "gray",
                                               size = 20)),
                              list(label = "C: Current age", shape = "icon",
                                   icon = list(code = "f007", color = "gray",
                                               size = 20))),
              useGroups = FALSE) %>%
    visHierarchicalLayout(sortMethod = "directed")

}
