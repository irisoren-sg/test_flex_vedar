require(data.table)
require(network)#for incidence matrix
require(intergraph) #for incidence matrix

combine_graphs <- function(g){
  #combine graphs in list of graphs g
  #if length of list == 1, return g[[1]]
  if(length(g) > 1){
  v <- map(g, 
           ~rbind(
                 igraph::as_data_frame(.x, 
                                       "vertices"
                                       )))
  #list to dataframe 
  # https://stackoverflow.com/questions/4227223/convert-a-list-to-a-data-frame
 v <- do.call(rbind.data.frame, v) %>% 
   remove_rownames() %>% 
   unique() 
  
  e <- map(g, 
           ~rbind(igraph::as_data_frame(.x, 
                                        "edges"))) 
  e <- do.call(rbind.data.frame, e) %>% 
    distinct()
  
  
  
  
  new_g <- igraph::graph_from_data_frame(
    e, 
    directed = TRUE, 
    vertices = v)
  }else{
    g[[1]]
  }
  
}
###########################################
extract_subgraph <- function(g, start_node, to = V(g),
                             direction = "out"){

    p <- igraph::all_simple_paths(g, from = start_node, 
                                  to = to, mode = direction)

  # get the vertices in the paths of interest
  vert <- unique(names(unlist(p)))
  
  #create the subgraph
  sub_g <- induced_subgraph(g, vert)
  sub_g
}

################################################################
extract_targets <- function(g, node_to_check, mode = "out", end = "head"){
  # extract the head/tail end vertices for node_to_check
  # return list of node names
  
  all_ends <- ends(g,
                   incident(g, node_to_check, mode = mode))
  if(end == "head"){
    all_ends[,2]
  }else if(end == "tail"){
    all_ends[,1]
  }else{
    stop("specify end as 'head' or 'tail'")
    
  }
}
#####################################

filter_graph_to_remove_repeats <- function(g, dat, sector_select, node_label){
  # nodes in g may have identical edges and be differentiated by only a digit
  # The function filters g to drop these repeated nodes 
  # If nodes have repeated names bit differing edges, all node versions are retained. This was done for simplicity
  # Inputs:
  #    g: graph
  #    dat:vd or vdt dat tibble
  #    sector_select: string. sector in dat that g belongs to
  #    node_label: column in dat used to label nodes in g
  # Output
  #    igraph object
  
  
  repeated_processes_check <- 
    #create tibble columns process_name (without digit suffix)|
    # data | heads_match (list) |tails_match  (list) |digit (list)
    check_edge_match_repeated_nodes(g, dat, sector_select) %>% 
    mutate(heads_match = unlist(heads_match), 
           tails_match = unlist(tails_match))
  
  processes_to_include <- repeated_processes_check %>% 
    filter(heads_match == TRUE &
             tails_match == TRUE) %>% 
    #create process column by pasting process_name and first entry of digit list
    mutate(process =  if_else(is.na(digit) == FALSE, 
                              paste(process_name, unlist(digit)[1], sep = ""), 
                              process_name)) %>% 
    select(process) 
  
  #handle the process with differing edges by unnesting data
  repeated_processes_with_diff_edges <- repeated_processes_check %>% 
    filter(heads_match == FALSE |
             tails_match == FALSE) 
  #if there is data, then this is appended to processes_to_include
  if(nrow(repeated_processes_with_diff_edges) > 0){
    repeated_processes_with_diff_edges <- 
      repeated_processes_with_diff_edges %>% 
      unlist_tibble_col(digit) %>% 
    #create process column of pasted process and digit
    mutate(process = map(process_name, ~unlist(.x)),
           process = paste(process_name, digit, sep = ""))
  
  processes_to_include <- processes_to_include %>% 
    bind_rows(repeated_processes_with_diff_edges %>% 
                select(process))
  }
  
  dat_unique <- dat %>% 
    filter(process %in% unique(processes_to_include$process))
  
  nodes_unique <- dat_unique %>% pull({{node_label}})
  
  vertices_to_drop <- V(g)$name[ 
    (V(g)$name %in% nodes_unique == F &
       (grepl("_process", V(g)$name) == F))]
  
  igraph::delete.vertices(g,  vertices_to_drop )  
}

###########################################
check_edge_match_repeated_nodes <- function(g, dat, sector_select){
  # procesess available in different periods may be distinguished by
  # digits
  # Check if processes with common string and different digits have 
  # the same connectivity in the graph
  # return tibble with logical for heads_match and tails_match
  # 
  # input:
  # g: graph of interest - assumes that V(g) labelled with process_description
  # dat: tibble of veda data from vdt
  # sector_select: sector in dat relevant to g
  
out <- dat %>% 
    # filter to data in g
    filter(process_description %in% V(g)$name) %>% 
    filter(sector == sector_select) %>% 
    # extract the digits from the process. process_name was not used as this was sometimes inconsistent across digit suffixes
    mutate(digit = str_extract(process, "[0-9]{2}")) %>% 
    #create columns with process name without the digit
    separate(process, sep = "([0-9]{2})",
             #column "other" specified in case anything 
             #follows digit
             into = c("process_name", "other"), 
             remove = F) %>% 
    #create a list of head nodes of all edges for each node
    mutate(heads = map(
      process_description, 
      ~extract_targets(g, .x, mode = "out", end = "head")
    ), 
    #create a list of head nodes of all edges for each node
    tails = map(
      process_description, 
      ~extract_targets(g, .x, mode = "in", end = "tail")
    ))  %>% 
    select(process_name, digit, heads, tails) %>% 
    #for each process_name (without the digit suffix), 
    #check if heads/tails are equal across suffixes
    #if only 1 process, return TRUE, if >2, return NA
    group_by(process_name) %>% 
    nest(data = c(digit, heads, tails)) %>% 
    # retain only unique rows
    mutate(data = map(data, ~unique(.x))) %>% 
    mutate(heads_match = map(data, 
           ~list_match_across_rows(.x, digit, heads)), 
           tails_match = map(data, 
                             ~list_match_across_rows(.x, digit,
                                                     tails)), 
           digit = map(data, 
                       ~unique(.x$digit))
    )
    
  out
  
}

list_match_across_rows <- function(dat, pivot_col, dat_col){
  #dat: tibble with pivot_col and dat_col
  #output: logical
  #TRUE if nrow(dat) == 1, 
  #if >2 rows check if all dat are identical
  
  
  new_cols <- dat %>% 
    select({{pivot_col}}) %>% 
    unique()
  
  if(nrow(dat) == 1){
    TRUE
  }else{
    #test if consecutive lists are identical. 
    # the returned logical will be a logical if all lists are identical
    # so we do not worry about comparing all lists, only consecutive lists
    t <- dat %>% 
      #lag data by one row
      mutate(shifted = as.vector(lag({{dat_col}}))) %>% 
      #lag introduces as null row
      filter(is.na(shifted) == F) %>% 
      # compare shifted data to dat_col
      mutate(identical_logical = map2({{dat_col}}, shifted, ~identical(.x, .y)))
    
    all(t$identical_logical == TRUE)  
    
  }
  
}

unlist_tibble_col <- function(dat, list_col){
  #unlist a tibble column to create new tibble rows with data repeated for non-listed entries
  
  # a tibble with list_col unlisted over rows by id
  
  unlisted <- dat %>%
               ungroup() %>%
               select({{list_col}}) 
  unlisted$id <-  seq.int(nrow(unlisted))
  unlisted <- unlisted %>% 
              unnest({{list_col}})


  #join unlisted to dat, and replace listed list_col with unnested list_col
  dat$id <- seq.int(nrow(dat))
  dat <- dat %>% 
    select(-{{list_col}}) %>% 
    left_join(unlisted, by = "id") %>% 
    select(-id)
  
  dat
}
####################

identify_terminal_nodes <- function(g){
  #return list of start and end nodes
  incidence_dat <- 
    #convert g to incidence matrix rows are the nodes
    as.matrix(asNetwork(g), matrix.type = "incidence") %>% 
    #transpose
    t() %>% 
    # and convert to tibble
    as_tibble(.name_repair = "minimal")
  
  #column names from vertices
  names(incidence_dat) = names(V(g))
  
  #count number of edges for each vertex
  edge_counts <- incidence_dat %>% 
    map(~sum(.x != 0)) %>% 
    as_tibble()
  #count incoming  edges for each vertex
  in_edge_counts <- incidence_dat %>% 
    map(~sum(.x == 1)) %>% 
    as_tibble()
  #count outgoing edges for each vertex
  out_edge_counts <- incidence_dat %>% 
    map(~sum(.x == -1)) %>% 
    as_tibble()
  
  edge_counts <- edge_counts %>% 
    bind_rows(out_edge_counts) %>% 
    bind_rows(in_edge_counts) 
  
  #start nodes have more than one out edge .x[2] and in edge=0 .x[3]
  start_nodes <- edge_counts %>% 
    map(function(.x).x[2]>=1 & .x[3]==0)
   #vice vera for end nodes     
  end_nodes <- edge_counts %>% 
    map(function(.x).x[2]==0 & .x[3]>=1)
  
  list( 
       start_nodes = names(start_nodes[which(start_nodes == T)]), 
       end_nodes = names(end_nodes[which(end_nodes == T)]))
        
}
####################

extract_subgraph_in_limited_graph <- function(g, node_select){
  #extract a subgraph for a node in a complex graph
  #
  # extract_subgraph() searches over entire graph
  # memory demanding for large graphs, so limit the search 
  #by searching from and to terminal nodes
  terminal_nodes <- identify_terminal_nodes(g)
  
  #graph for all edges into node_select
  in_graph <- map(terminal_nodes$start_nodes,~extract_subgraph(g, 
                                                                     start_node = node_select,
                                                                     to = .x,
                                                                     direction = "in"))
  #graph for all edges out of node_select
  out_graph <- map(terminal_nodes$end_nodes,~extract_subgraph(g, 
                                                               start_node = node_select,
                                                               to = .x,
                                                               direction = "out"))
  #combine lists
  all_g <- c(in_graph, out_graph) %>% 
    # into graph
    combine_graphs() 
  
}

