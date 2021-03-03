#'---
#'title: Lineage Naming Function
#'author: Kathryn Campbell
#'date: 03/03/2021
#'---
#'
#'
#'**LINEAGE NAMING**
#'=================
#'
#'File: lineage_naming.r
#'=====================

lineage_naming <- function(sequence_data, metadata, node_data, tree) {
  
sequence_data$previous <- NA
for (i in 1:length(sequence_data$ID)) {
  sequence_data$previous[i]<-
    metadata$alignment.displayName[which(metadata$ID == sequence_data$ID[i])]
}

previous_assignments<-data.frame(assignment = unique(sequence_data$previous), node = NA)

node_data$previous<-NA

for (i in 1:length(node_data$Node)) {
  clades<-unique(sequence_data$previous[
    which(sequence_data$ID %in% tree$tip.label[c(unlist(
      Descendants(tree, node_data$Node[i], type = "tips")))])])
  
  node_data$previous[i]<-
    paste(c(clades), collapse = ", ")
  
}

for (i in 1:length(previous_assignments$assignment)) {
  previous_assignments$node[i]<-which(node_data$previous == previous_assignments$assignment[i])[1]
  previous_assignments$assignment[i]<-str_split(previous_assignments$assignment, " ")[[i]][2]
}

previous_assignments<-previous_assignments[-c(which((is.na(previous_assignments$node) | is.na(previous_assignments$assignment)))),]
possible_names<-data.frame(names = rep(previous_assignments$assignment, 9))
previous_assignments$assignment<-paste(previous_assignments$assignment, "_A1", sep = "")

for (i in 1:length(previous_assignments$assignment)) {
  node_data$cluster[previous_assignments$node[i]]<-previous_assignments$assignment[i]
}

node_data$cluster[1]<-"A1"
node_data$test <- NA
problem_names<-data.frame(letters = c("A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "I1"))
possible_names<-possible_names[order(possible_names$names),]
possible_names<-paste(possible_names, problem_names$letters, sep = "_")

for (i in 1:length(node_data$Node)) {
  test<-which(node_data$Node %in% descendants(tree, node_data$Node[i], type = "all", ignore.tip = T))
  node_data$test[c(test)] <- paste(node_data$cluster[i], ".1", sep = "")
  node_data$test<-str_replace(node_data$test, "A1\\..\\..\\..", "B1")
  node_data$test<-str_replace(node_data$test, "B1\\..\\..\\..", "C1")
  node_data$test<-str_replace(node_data$test, "C1\\..\\..\\..", "D1")
  node_data$test<-str_replace(node_data$test, "D1\\..\\..\\..", "E1")
  node_data$test<-str_replace(node_data$test, "E1\\..\\..\\..", "F1")
  node_data$test<-str_replace(node_data$test, "F1\\..\\..\\..", "G1")
  node_data$test<-str_replace(node_data$test, "G1\\..\\..\\..", "H1")
  
  majors<-which(grepl("_", node_data$test))
  node_data$cluster[c(majors)] <- node_data$test[c(majors)]
  
  for (k in 1:length(possible_names)) {
    if (length(which(node_data$cluster == possible_names[k]))>1) {
      problems<-which(node_data$cluster == possible_names[k])
      problems<-problems[-c(1)]
      y=1
      for (a in 1:length(problems)) {
        letter<-which(problem_names$letters == (str_split(node_data$cluster[problems[a]], "_")[[1]][2]))
        node_data$cluster[problems[a]]<-paste((str_split(node_data$cluster[problems[a]], "_")[[1]][1]), problem_names$letters[(letter+y)], sep = "_")
        y = y+1
      }
    }
  }
  duplicates<-unique(node_data$cluster[duplicated(node_data$cluster)])
  problems<-duplicates[which(str_count(duplicates, pattern = "\\.") == 0)]
  duplicates<-duplicates[which(str_count(duplicates, pattern = "\\.") != 0)]
  
  for (i in 1:length(duplicates)) {
    test<-which(node_data$cluster == duplicates[i])
    test<-test[-c(1)]
    x<-1
    for (j in 1:length(test)) {
      name<-unlist(str_split(node_data$cluster[test[j]], "\\."))
      name[length(name)]<-x+as.integer(name[length(name)])
      x<-(x+1)
      node_data$cluster[test[j]]<-paste(c(name), collapse='.' )
    }
  }
}

unclassified<-which(!grepl("_", node_data$cluster))
unclassified<-unclassified[c(-1)]
for (i in 1:length(node_data$Node)) {
  test<-which(node_data$Node %in% descendants(tree, node_data$Node[i], type = "all", ignore.tip = T))
  node_data$test[c(test)] <- paste(node_data$cluster[i], ".1", sep = "")
  node_data$test<-str_replace(node_data$test, "A1\\..\\..\\..", "B1")
  node_data$test<-str_replace(node_data$test, "B1\\..\\..\\..", "C1")
  node_data$test<-str_replace(node_data$test, "C1\\..\\..\\..", "D1")
  node_data$test<-str_replace(node_data$test, "D1\\..\\..\\..", "E1")
  node_data$test<-str_replace(node_data$test, "E1\\..\\..\\..", "F1")
  node_data$test<-str_replace(node_data$test, "F1\\..\\..\\..", "G1")
  node_data$test<-str_replace(node_data$test, "G1\\..\\..\\..", "H1")
  
  node_data$cluster[unclassified]<-node_data$test[unclassified]
  
  for (v in 1:length(problem_names$letters)) {
    if (length(which(node_data$cluster == problem_names$letters[v]))>1) {
      problems<-which(node_data$cluster == problem_names$letters[v])
      problems<-problems[-c(1)]
      y=1
      for (f in 1:length(problems)) {
        letter<-which(problem_names$letters == (node_data$cluster[problems[f]]))
        node_data$cluster[problems[f]]<-problem_names$letters[(letter+y)]
        y = y+1
      }
    }
  }
  duplicates<-unique(node_data$cluster[duplicated(node_data$cluster)])
  problems<-duplicates[which(str_count(duplicates, pattern = "\\.") == 0)]
  duplicates<-duplicates[which(str_count(duplicates, pattern = "\\.") != 0)]
  
  for (i in 1:length(duplicates)) {
    test<-which(node_data$cluster == duplicates[i])
    test<-test[-c(1)]
    x<-1
    for (j in 1:length(test)) {
      name<-unlist(str_split(node_data$cluster[test[j]], "\\."))
      name[length(name)]<-x+as.integer(name[length(name)])
      x<-(x+1)
      node_data$cluster[test[j]]<-paste(c(name), collapse='.' )
    }
  }
}

node_data<-node_data[,-c(6,7)]
return(node_data)
}
