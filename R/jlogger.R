#logger that uses cat instead of print
#depends on util.R
if(!"CFG" %in% ls()) CFG = NULL
JLogger = setRefClass("JLogger", fields = list(m_level = "integer", m_file = "character"))
JLOGGER.TRACE = 1L
JLOGGER.DEBUG = 2L
JLOGGER.INFO = 3L
JLOGGER.WARN = 4L
JLOGGER.ERROR = 5L
JLOGGER.FATAL = 6L
JLOGGER.LEVELS = c("TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL")
JLOGGER.DEFAULT_LEVEL = JLOGGER.DEBUG
JLoggerFactory <- function(name, ..., reset = FALSE)
{
    if(!"JLOGGER.ENV" %in% ls(.GlobalEnv)) JLOGGER.ENV <<- new.env()
    if(!name %in% ls(JLOGGER.ENV) || reset) assign(name, JLogger(name, ...), JLOGGER.ENV)
    get(name, JLOGGER.ENV)
}

JLoggerReset <- function()
{
    rm(list = ls(env = JLOGGER.ENV), envir = JLOGGER.ENV)
}

SetJLogger <- function(name, jlogger, ...)
{
    if(!"JLOGGER.ENV" %in% ls(.GlobalEnv)) JLOGGER.ENV <<- new.env()
    assign(name, jlogger, JLOGGER.ENV)
}

#We want the JLogger to be a ref class so the prefix can be changed without having to be propagated
#m_files can be several file if we want to write to several connections at the same time
JLogger <- setRefClass("JLogger",
                       fields = list(m_files = "character", m_level = "integer", m_name = "character", m_prefix = "character"))

JLOGGER.init <- function(name = "",
                         files = "",
                         prefix = name,
                         level,
                         logconfig = CFG$jlogger)
{
    m_name <<- name
    m_files <<- files
    if(missing(level))level = JLOGGER.getlevel(name, logconfig)
    m_level <<- level
    m_prefix <<- prefix
}

JLogger$methods(initialize = JLOGGER.init)

#We may want to flush
JLOGGER.flush <- function(jlfile)
{
   if(jlfile != "") cat("", file = jlfile)
}

JLOGGER.getlevel <- function(name, logconfig)
{
    if(missing(logconfig) || is.null(logconfig)) return(JLOGGER.DEFAULT_LEVEL)
    level = 0L
    for(lv in JLOGGER.LEVELS)
    {
        level = level + 1L
        lvnode = logconfig[lv]
        if(is.null(lvnode)) next
        if(mgrep(lvnode, grep_fun = grepl, name)) return(level)        
    }
    return(level + 1)#won't log anything
}

#Not a true set will return a copy. Can override a few value but not the name and file
jlset <- function(jlogger, level, ...)
{
    if(missing(level))level=jlogger@level
    return(JLogger(jlogger$m_name, jlogger$m_file, level = level, ...))
}

#Is the logger quiet for the given level
JLOGGER.jlquiet <- function(jlogger, level, ...)
{
    if(is.null(jlogger)) return(TRUE)
    level < jlogger$m_level
}

#To use for object that can't be printed with cat. Uses write.table internally
JLOGGER.jlwrite <- function(jlfile, level, prefix, data, ..., endline = "\n")
{
    cat(as.character(Sys.time()), JLOGGER.LEVELS[level], prefix, ":", endline, file = jlfile, append = TRUE)
    suppressWarnings(write.table(data, ..., file = jlfile, append = TRUE))#Warns about appending column names to a file
    cat(endline, file = jlfile, append = TRUE)
}

#To print complicated objects not handle by cat to the console
JLOGGER.jlprint <- function(jlfile, level, prefix, data,..., endline = "\n")
{
    #This prints to the console so only id jlfile == ""
    if(jlfile != "") return()
    cat(as.character(Sys.time()), JLOGGER.LEVELS[level], prefix, ":", endline, file = jlfile, append = TRUE)
    print(data)
    cat(endline, file = jlfile, append = TRUE)
}

JLOGGER.jlog <- function(jlfile,
                         level,
                         ...,
                         prefix,
                         prechar = "",
                         endline = "\n",
                         cat_fun = cat)
{
    if(prechar != "") cat(prechar, file = jlfile, append = TRUE)
    cat_fun(as.character(Sys.time()), JLOGGER.LEVELS[level], prefix, ":", ..., endline, file = jlfile, append = TRUE)
}

JLOGGER.do <- function(jlogger, level, log_fun, ..., prefix = jlogger$m_prefix )
{
    if(JLOGGER.jlquiet(jlogger, level, ...)) return()
    lapply(jlogger$m_files, log_fun, prefix = prefix, level = level, ...)
    invisible()
}

jl.log.trace <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.TRACE, JLOGGER.jlog, ...)
jl.log.debug <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.DEBUG, JLOGGER.jlog, ...)
jl.log.info <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.INFO, JLOGGER.jlog, ...)
jl.log.warn <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.WARN, JLOGGER.jlog, ...)
jl.log.error <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.ERROR, JLOGGER.jlog, ...)
jl.log.fatal <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.FATAL, JLOGGER.jlog, ...)

jl.write.trace <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.TRACE, JLOGGER.jlwrite, ...)
jl.write.debug <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.DEBUG, JLOGGER.jlwrite, ...)
jl.write.info <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.INFO, JLOGGER.jlwrite, ...)
jl.write.warn <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.WARN, JLOGGER.jlwrite, ...)
jl.write.error <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.ERROR, JLOGGER.jlwrite, ...)
jl.write.fatal <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.FATAL, JLOGGER.jlwrite, ...)

jl.print.trace <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.TRACE, JLOGGER.jlprint, ...)
jl.print.debug <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.DEBUG, JLOGGER.jlprint, ...)
jl.print.info <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.INFO, JLOGGER.jlprint, ...)
jl.print.warn <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.WARN, JLOGGER.jlprint, ...)
jl.print.error <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.ERROR, JLOGGER.jlprint, ...)
jl.print.fatal <- function(jlogger, ...) JLOGGER.do(jlogger, JLOGGER.FATAL, JLOGGER.jlprint, ...)


setGeneric("jlog.trace", function(jlogger, ...) standardGeneric("jlog.trace"))
setMethod("jlog.trace", "JLogger", jl.log.trace)
setGeneric("jlwrite.trace", function(jlogger, ...) standardGeneric("jlwrite.trace"))
setMethod("jlwrite.trace", "JLogger", jl.write.trace)
setGeneric("jlprint.trace", function(jlogger, ...) standardGeneric("jlprint.trace"))
setMethod("jlprint.trace", "JLogger", jl.print.trace)

setGeneric("jlog.debug", function(jlogger, ...) standardGeneric("jlog.debug"))
setMethod("jlog.debug", "JLogger", jl.log.debug)
setGeneric("jlwrite.debug", function(jlogger, ...) standardGeneric("jlwrite.debug"))
setMethod("jlwrite.debug", "JLogger", jl.write.debug)
setGeneric("jlprint.debug", function(jlogger, ...) standardGeneric("jlprint.debug"))
setMethod("jlprint.debug", "JLogger", jl.print.debug)

setGeneric("jlog.info", function(jlogger, ...) standardGeneric("jlog.info"))
setMethod("jlog.info", "JLogger", jl.log.info)
setGeneric("jlwrite.info", function(jlogger, ...) standardGeneric("jlwrite.info"))
setMethod("jlwrite.info", "JLogger", jl.write.info)
setGeneric("jlprint.info", function(jlogger, ...) standardGeneric("jlprint.info"))
setMethod("jlprint.info", "JLogger", jl.print.info)

setGeneric("jlog.warn", function(jlogger, ...) standardGeneric("jlog.warn"))
setMethod("jlog.warn", "JLogger", jl.log.warn)
setGeneric("jlwrite.warn", function(jlogger, ...) standardGeneric("jlwrite.warn"))
setMethod("jlwrite.warn", "JLogger", jl.write.warn)
setGeneric("jlprint.warn", function(jlogger, ...) standardGeneric("jlprint.warn"))
setMethod("jlprint.warn", "JLogger", jl.print.warn)

setGeneric("jlog.error", function(jlogger, ...) standardGeneric("jlog.error"))
setMethod("jlog.error", "JLogger", jl.log.error)
setGeneric("jlwrite.error", function(jlogger, ...) standardGeneric("jlwrite.error"))
setMethod("jlwrite.error", "JLogger", jl.write.error)
setGeneric("jlprint.error", function(jlogger, ...) standardGeneric("jlprint.error"))
setMethod("jlprint.error", "JLogger", jl.print.error)

setGeneric("jlog.fatal", function(jlogger, ...) standardGeneric("jlog.fatal"))
setMethod("jlog.fatal", "JLogger", jl.log.fatal)
setGeneric("jlwrite.fatal", function(jlogger, ...) standardGeneric("jlwrite.fatal"))
setMethod("jlwrite.fatal", "JLogger", jl.write.fatal)
setGeneric("jlprint.fatal", function(jlogger, ...) standardGeneric("jlprint.fatal"))
setMethod("jlprint.fatal", "JLogger", jl.print.fatal)

jl.flush <- function(jlogger)
{
    lapply(jlogger$m_files, JLOGGER.flush)
    invisible()
}

setGeneric("jlflush", function(jlogger) standardGeneric("jlflush"))
setMethod("jlflush", "JLogger", jl.flush)

ljlog.trace <- function(jlogger, ...) lapply(jlogger, jlog.trace, ...)
ljlwrite.trace <- function(jlogger, ...) lapply(jlogger, jlwrite.trace, ...)
ljlprint.trace <- function(jlogger, ...) lapply(jlogger, jlprint.trace, ...)

ljlog.debug <- function(jlogger, ...) lapply(jlogger, jlog.debug, ...)
ljlwrite.debug <- function(jlogger, ...) lapply(jlogger, jlwrite.debug, ...)
ljlprint.debug <- function(jlogger, ...) lapply(jlogger, jlprint.debug, ...)

ljlog.info <- function(jlogger, ...) lapply(jlogger, jlog.info, ...)
ljlwrite.info <- function(jlogger, ...) lapply(jlogger, jlwrite.info, ...)
ljlprint.info <- function(jlogger, ...) lapply(jlogger, jlprint.info, ...)

ljlog.warn <- function(jlogger, ...) lapply(jlogger, jlog.warn, ...)
ljlwrite.warn <- function(jlogger, ...) lapply(jlogger, jlwrite.warn, ...)
ljlprint.warn <- function(jlogger, ...) lapply(jlogger, jlprint.warn, ...)

ljlog.error <- function(jlogger, ...) lapply(jlogger, jlog.error, ...)
ljlwrite.error <- function(jlogger, ...) lapply(jlogger, jlwrite.error, ...)
ljlprint.error <- function(jlogger, ...) lapply(jlogger, jlprint.error, ...)

ljlog.fatal <- function(jlogger, ...) lapply(jlogger, jlog.fatal, ...)
ljlwrite.fatal <- function(jlogger, ...) lapply(jlogger, jlwrite.fatal, ...)
ljlprint.fatal <- function(jlogger, ...) lapply(jlogger, jlprint.fatal, ...)

setMethod("jlog.trace", "list", ljlog.trace)
setMethod("jlwrite.trace", "list", ljlwrite.trace)
setMethod("jlprint.trace", "list", ljlprint.trace)

setMethod("jlog.debug", "list", ljlog.debug)
setMethod("jlwrite.debug", "list", ljlwrite.debug)
setMethod("jlprint.debug", "list", ljlprint.debug)

setMethod("jlog.info", "list", ljlog.info)
setMethod("jlwrite.info", "list", ljlwrite.info)
setMethod("jlprint.info", "list", ljlprint.info)

setMethod("jlog.warn", "list", ljlog.warn)
setMethod("jlwrite.warn", "list", ljlwrite.warn)
setMethod("jlprint.warn", "list", ljlprint.warn)

setMethod("jlog.error", "list", ljlog.error)
setMethod("jlwrite.error", "list", ljlwrite.error)
setMethod("jlprint.error", "list", ljlprint.error)

setMethod("jlog.fatal", "list", ljlog.fatal)
setMethod("jlwrite.fatal", "list", ljlwrite.fatal)
setMethod("jlprint.fatal", "list", ljlprint.fatal)

ljlflush <- function(jlogger) lapply(jlogger, jlflush)
setMethod("jlflush", "list", ljlflush)

set_log_files <- function(logger,
                          files)
{
    if(is.character(logger)) logger = JLoggerFactory(logger)
    logger$m_files = unique(files)
}

add_log_files <- function(logger,
                          files)
{
   if(is.character(logger)) logger = JLoggerFactory(logger)
   set_log_files(logger, c(logger$m_files, files))
}

###### utility function for logging more information ####

## Returns the name of the function it was called from
fname <- function(offset = 1)
{
    as.character(sys.calls()[[sys.nframe() - offset]])[1]    
}

jlfname <- function() fname(8)
