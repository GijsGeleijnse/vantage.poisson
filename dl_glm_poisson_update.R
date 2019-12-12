library(mltools)
library(rjson)
library(dplyr)
library(stats4)

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


pre_process_data <- function(df, expl_vars, response_col) {
  # Split dataframe into explanatory and response variables,  add intercept
  X <- df[, expl_vars]
  Y <- df[,response_col]

  return(list(
    X = X,
    Y = Y
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

  return(-s)
}


dl_compute_loglikelihood <- function(beta, client, call.method) {
  aggregates <- call.method(client, "perform_iteration", expl_vars, response_col, beta)

  s <- 0
  for (k in 1:length(aggregates)){
    s <- s + aggregates[[k]]
  }

  return(s)
}

dl_glm_poisson <- function(client, expl_vars, response_col, call.method) {

    # Define p0 as a closure (with access to client, expl_vars, ...)
    p0 <- function() {
        vars <- as.list(environment())
        do.call(
            dl_compute_loglikelihood,
            list(
                beta=vars,
                client=client,
                call.method=call.method
            )
        )
    }

    # Create a list with keys beta*=0.
    initial <- create.betas(length(expl_vars))

    # Override the function definition of p0, explicitly setting
    # the arguments to be equal to the list `initial`
    formals(p0) <- initial

    res <- stats4::mle(p0, start=initial)
    return(res)
}


dl_poisson.mock <- function(df, expl_vars, response_col, splits=5) {

  datasets <- list()

  for (k in 1:splits) {
    datasets[[k]] <- df[seq(k, nrow(df), by=splits), ]
  }

  # No more global assignment ...
  client <- MockClient(datasets)
  results <- dl_glm_poisson(client, expl_vars, response_col, mock.call)

  return(results)
}

mock.example <- function(){
    library(data.table)

    # Load table and replace column "prog" with a factor
    p <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
    p <- within(p, {
      prog <- factor(
          prog,
          levels=1:3,
          labels=c("General", "Academic", "Vocational")
      )
    })

    pp <- one_hot(as.data.table(p))
    pp <- as.data.frame(pp)
    pp$intercept <- rep(1,nrow(pp))

    vv <- dl_poisson.mock(
        df=pp,
        expl_vars=c("prog_Academic","prog_Vocational","math","intercept"),
        response_col = "num_awards"
    )

  summary(m1 <- glm(num_awards ~ prog_Academic + prog_Vocational + math, family="poisson", data=pp))
}

create.betas <- function(nr) {
    my_list <- list()
    for(i in 1:nr){
        varname <- paste("beta", i-1, sep='')
        my_list[[varname]] <- 0

    }

    return(my_list)
}

# Run the example when source-ing ...
mock.example()
