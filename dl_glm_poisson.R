#

#' ----------------------------------------------------------------------------
#' description:
#'   Implementation of the distributed Poisson Regression 
#'   This implementation can be used with VANTAGE
#'   
#'   
#' author:
#'   Melle Sieswerda <m.sieswerda@iknl.nl>
#'   Gijs Geleijnse <g.geleijnse@iknl.nl>
#' date: November 29, 2019
#' license: MIT License
#' ----------------------------------------------------------------------------




library(mltools)

library(rjson)
library(dplyr)



########### FUNCTIONS copied from dl_coxph
############################################


#' Entrypoint when excecuting this script using Rscript
#'
#' Wraps the docker input/output for `dispatch_RPC()`.
#' Deserialization/serialization is performed in `dipatch_RPC()` to enable
#' testing.
docker.wrapper <- function() {
  database_uri <- Sys.getenv("DATABASE_URI")
  writeln(sprintf("Using '%s' as database", database_uri))
  df <- read.csv(database_uri)
  
  # Read the contents of file input.txt into 'input_data'
  writeln("Loading input.txt")
  filename <- 'input.txt'
  input_data <- readChar(filename, file.info(filename)$size)
  
  writeln("Dispatching ...")
  result <- dispatch_RPC(df, input_data)
  
  # Write result to disk
  writeln("Writing result to disk .. ")
  writeLines(result, "output.txt")
  
  writeln("")
  writeln("[DONE!]")
}


# ******************************************************************************
# ---- main() ----
# ******************************************************************************
if (!interactive()) {
  docker.wrapper()
}


#' Mock an RPC call to all sites.
#'
#' Params:
#'   client: ptmclient::Client instance.
#'   method: name of the method to call on the distributed learning
#'           infrastructure
#'   ...: (keyword) arguments to provide to method. The arguments are serialized
#'        using `saveRDS()` by `create_task_input()`.
#'
#' Return:
#'   return value of called method
mock.call <- function(client, method, ...) {
  
  writeln(sprintf('** Mocking call to "%s" **', method))
  datasets <- client$datasets
  input_data <- create_task_input(method, ...)
  input_data <- toJSON(input_data)
  
  # Create a list to store the responses from the individual sites
  results <- list()
  
  # Mock calling the RPC method on each site
  for (k in 1:length(datasets)) {
    result <- dispatch_RPC(datasets[[k]], input_data)
    results[[k]] <- readRDS(textConnection(result))
  }
  
  writeln()
  return(results)
}



#' Create a data structure used as input for a call to the distributed
#' learning infrastructure.
create_task_input = function(method, ...) {
  # Construct the input_data list from the ellipsis.
  arguments <- list(...)
  
  if (is.null(names(arguments))) {
    args <- arguments
    kwargs <- list()
    
  } else {
    args <- arguments[names(arguments) == ""]
    kwargs <- arguments[names(arguments) != ""]
  }
  
  # Serialize the argument values to ASCII
  fp <- textConnection("arg_data", open="w")
  saveRDS(args, fp, ascii=T)
  close(fp)
  
  # Serialize the keyword argument values to ASCII
  fp <- textConnection("kwarg_data", open="w")
  saveRDS(kwargs, fp, ascii=T)
  close(fp)
  
  # Create the data structure
  input_data <- list(
    method=method,
    args=arg_data,
    kwargs=kwarg_data
  )
  
  return(input_data)
}


# ******************************************************************************
# ---- Infrastructure functions ----
# ******************************************************************************

#' Run the method requested by the server
#'
#' Params:
#'   df: data frame containing the *local* dataset
#'   input_data: string containing serialized JSON; JSON should contain
#'               the keys 'method', 'args' and 'kwargs'
#'
#' Return:
#'   Requested method's output
dispatch_RPC <- function(df, input_data) {
  # Determine which method was requested and combine arguments and keyword
  # arguments in a single variable
  input_data <- fromJSON(input_data)
  method <- sprintf("RPC_%s", input_data$method)
  
  input_data$args <- readRDS(textConnection(input_data$args))
  input_data$kwargs <- readRDS(textConnection(input_data$kwargs))
  
  args <- c(list(df), input_data$args, input_data$kwargs)
  
  # Call the method
  writeln(sprintf("Calling %s", method))
  result <- do.call(method, args)
  
  # Serialize the result
  writeln("Serializing result")
  fp <- textConnection("result_data", open="w")
  saveRDS(result, fp, ascii=T)
  close(fp)
  result <- result_data
  
  return(result)
}


#####################################################################
##################### Poisson Regression
#####################################################################




#### helper functions at the node ############

pre_process_data <- function(df, expl_vars, response_col) {
 
  
  # Split dataframe into explanatory and response variables,  add intercept
  X <- df[, expl_vars]
 
  Y <- df[,response_col]
  
  
  
  return(list(
    X = X, Y = Y
  ))
}







# ******************************************************************************
# ---- RPC entry points -- node ----
# ******************************************************************************

RPC_perform_iteration <- function(df, expl_vars, response_col, beta) {
  
  data <- pre_process_data(df, expl_vars, response_col)
  
  X <- data$X
  Y <- data$Y
  
  Y <- as.integer(Y)
  
  ## Sum_i_in_m  (y_i beta x_i - exp(beta * x_i) - log(Y_i!) )
  s <- 0
  for(i in 1:nrow(X)){
    t <- sum(X[i,] * beta)
    u <- Y[i] * t
    u <- u - exp(t)
    u <- u - log(factorial(Y[i]))
    s <- s + u
    
    
    
  }
 
  -s
}


dl_compute_loglikelihood <- function(beta){
  
  s <- 0
  aggregates <- call.method(client, "perform_iteration", expl_vars, response_col, beta)
  for (k in 1:length(aggregates)){
    s <- s + aggregates[[k]]
    
  }
  s
}

p1 <-  function(beta0){
  
  dl_compute_loglikelihood(beta=c(beta0))
  
}


p2 <-  function(beta0,beta1){
  
  dl_compute_loglikelihood(beta=c(beta0,beta1))
  
}

p3 <-  function(beta0,beta1,beta2){
  
  dl_compute_loglikelihood(beta=c(beta0,beta1,beta2))
  
}



p4 <-  function(beta0,beta1,beta2,beta3){
 
  dl_compute_loglikelihood(beta=c(beta0,beta1,beta2,beta3))
  
}

p5 <-  function(beta0,beta1,beta2,beta3,beta4){
  
  dl_compute_loglikelihood(beta=c(beta0,beta1,beta2,beta3,beta4))
  
}

p6 <-  function(beta0,beta1,beta2,beta3,beta4,beta5){
  
  dl_compute_loglikelihood(beta=c(beta0,beta1,beta2,beta3,beta4,beta5))
  
}



dl_glm_poisson <- function(){
  ## assumption: intercept is  one of the expl_vars
  
  m <- length(expl_vars)
  
  if(m==1){ 
    
    res <- mle(p1,start=list(beta0=0))
  } else
  if(m==2){ 
    
    res <- mle(p2,start=list(beta0=0,beta1=0))
  } else
  
  if(m==3){ 
    
    res <- mle(p3,start=list(beta0=0,beta1=0,beta2=0))
  } else
  if(m==4){ 
  
  res <- mle(p4,start=list(beta0=0,beta1=0,beta2=0,beta3=0))
  } else
  if(m==5){ 
      
      res <- mle(p5,start=list(beta0=0,beta1=0,beta2=0,beta3=0, beta4=0))
    }
  
  
  return(res)
}








dl_poisson.mock <- function(df, expl_vars, response_col, splits=5) {
  
  datasets <- list()
  
  for (k in 1:splits) {
    datasets[[k]] <- df[seq(k, nrow(df), by=splits), ]
  }
  
  client <<- MockClient(datasets)
  expl_vars <<- expl_vars
  response_col <<- response_col
  call.method<<- mock.call
  
  results <- dl_glm_poisson()
  return(results)
}

mock.example <- function(){
  library(data.table)
  p <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
  
  p <- within(p, {
    prog <- factor(prog, levels=1:3, labels=c("General", "Academic", 
                                              "Vocational"))
    
  })
  
  
  
  with(p, tapply(num_awards, prog, function(x) {
    sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
  }))
  
  pp <- one_hot(as.data.table(p))
  pp <- as.data.frame(pp)
  pp$intercept <- rep(1,nrow(pp))
  
  dl_poisson.mock(df=pp,expl_vars=c("prog_Academic","prog_Vocational","intercept"),response_col = "num_awards")
  
}
