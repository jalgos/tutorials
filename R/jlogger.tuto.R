## Gracefully handling logging is key when developping big projects
## Logging has several purposes:
##   - outputing the result of a computation
##   - printing warning when the program does something the user does not necessarily expect
##   - printing detailed explanations of a program's error/ crash
##   - Providing the user with the current state of computation (e.g. "step i over n in function f"
##   - Helping debugging by tracing the successive calls and variable values to identify a bug

## Furthermore we may want to output the logs either directly to the console or write them in files to be stored.

## R has a few functions to handle printing objects, namely `print`, `message` and `cat`. However they have limited functionality.
## That is why we built the 'JLogger' class.
## JLogger is built on top of cat and allows one to output to several files or to the console easily. JLogger statements can be turn on and off as it implements the well known concept of logging level (TRACE, DEBUG, INFO, WARNING, ERROR, FATAL)

## JLogger are R5 objects. This is one of the few occurences of R5 objects as their use is strongly discouraged.
## Mutable objects where chosen here to persist changes made locallyto a logger. For example if a condition is met in a function the logging level or logging file of a logger can be changed affecting the rest of the program.

##################################################
## I) Use of logging
##################################################

## Logging should be used to give the user key information about the current of state of computation.
## It is really important to print information about the size of the inputs, the current step in a recursive algorithm, measures about the outputs, edge cases for example.

## Logging should be implemented from the start of the project.


## It is very common, especially early in the development of a project, that we need to trace a computation in its tiniest details. This situation arises for debugging purposes mostly.

## The way one will go about it would be to put 'print'/ 'cat' statements in points of interest.

## Let say we want to debug the following function we can insert a few cat statements
f <- function()
{
    X = rnorm(1000)
    cat("Maximum value of X is:", max(X), "\n")
    D = rep(0, length(X) + 1)
    for(i in 1:length(X))
    {
        cat("Value of X at step i:", X[i], "\n")
        D[i + 1] = cos(D[i] + X[i])
    }
    cat("Median of D is:", median(D), "\n")
    sum(D)
}

f()

## We see have a lot of information about the computation and monitor it easily now.
## The problem is when we fixed the function the call in the for loop is of no interest but yet will take most of the logging space, making it hard to find relevant information for further debugging.

## A solution would be to remove or comment the call to 'cat' in the loop and proceed with the new version.
## However if in the future our function gives us an unexpected output that requires an in-depth look at the function's computation wewould have to put the statement back.
## It can be easy to handle for a small number of functions but becomes a nightmare as the project gets bigger.

## JLogger is a perfect fit here. Call to the logger can be turned on and off according to the logging level of the logger

##################################################
## II) Using JLogger
##################################################

## First you need to source the file "jlogger.R"
## The file can be found in the git project 'util' along with other utility R functions. url: 'git@datasaiyan.com:/home/git/repos/utils/util.R.git'
## Jlogger requires you to source the file "util.R" as well (uses function 'mgrep')

## JLogger objects are created/ retrieved through the JLoggerFactory
## Each logger is identified by a string id:

logger = JLoggerFactory("jalgos")

## If the logger named 'jalgos' was not created before the above call will create a new instance. Otherwise it will return the object previously created
## This makes it easy to use the logger from within any function without having to pass it around.

## JLoggers can filter which statements they output to the log file (or console). It is achieved through the logging level
## 6 levels are implemented
JLOGGER.LEVELS

## The above vectors defines the logging priority. For example if a logger level is set to 'JLOGGER.TRACE' any message will be output to the log file.
## Conversely if it is set to 'JLOGGER.FATAL' only "fatal" messages will be print out.
## Let's illustrate this behavior:

logger$m_level = JLOGGER.INFO

## the syntax for logging a message is:
## jlog.$priority($logger, $object1, $object2, ...., $objectn)
## The logger will check if the call has a priority greater than its priority and output a message that looks like "date time loglevel prefix: $object1, $object2, ...., $objectn \n"
## No need to give the eol character or separate objects with spaces

jlog.info(logger, "This is a test")
## The logger gives a bit more context that can be usefull in debugging or monitoring execution.

## Since the logging level is set to info the following message won't be logged

jlog.trace(logger, "Tracing")



my_heavy_comp <- function()
{
    logger = JLoggerFactory("jalgos")
    jlog.info(logger, "Function:", jlfname(), "This computation is gonna take a long time")
    X = rnorm(1000)
    jlog.debug(logger, "Max of X:", max(X))
    D = 0
    for(x in X)
    {
        jlog.trace(logger, "I want to know exactly what is going on at each step of the loop:", "x:", x, "D:", D)
        D = sin(x + D)
        if(abs(D) < .01) jlog.warn(logger, "Function:", jlfname(), "D is small:", D, ", you should know")
    }
    jlog.info(logger, "Function:", jlfname(), "Computation is over. Return value:", D)
    D
}

## setting logger level to info
logger = JLoggerFactory("jalgos")
logger$m_level = JLOGGER.INFO
my_heavy_comp()

## Now to warning:
logger$m_level = JLOGGER.WARN
my_heavy_comp()

## And now to trace:

logger$m_level = JLOGGER.TRACE
my_heavy_comp()

## The function 'jlfname' is specifically designed to print the function name from within a jlog.* call. It works only from there. Trying to use it elsewhere in the code will give unpredictable results.
## 'jlfname' calls 'fname' that relies on the current state of the call stack. Because of promissed evaluation the stack will change at the time 'fname' is evaluated.
## 'fname' called from the body of a function gives the function name if evaluated in the body.

fblah <- function() fname()
fblah()

## This example illustrates how messaging can be turned on and off easily
## One should not refrain from using the logger. Appropriately logging relevant information really improves post-mortem debugging and greatly speeds up development

## We provide two other ways of logging objects: 'jlprint.*' and 'jlwrite.*'
## 'jlprint' is a wrapper around 'print' and 'jlwrite' is a wrapper around 'write.table'. This allows you to log more complex objects that are not handled by cat.
m = matrix(rnorm(100), 10, 10)
jlog.info(logger, m) ## Fails
jlprint.info(logger, m)
d = data.frame(i = 1:10, x = rnorm(10))
jlwrite.info(logger, d)

##################################################
## III) Configuring JLogger
##################################################

## Logging levels 
## Logging levels can be set directly:
logger = JLoggerFactory("jalgos")
logger$m_level = JLOGGER.INFO

## It can also be set through the function 'set_logging_level'
set_logging_level(logger, JLOGGER.INFO)

## 'set_logging_level' handles logger objects and logger names
set_logging_level("jalgos", JLOGGER.INFO)

## Output files
## The logger can output its messages to one or several files
## The logging files can be access and overriden through the member 'm_files'

logger$m_files
## This field is a character vector of the log file names.
## The only exception is that the filename "" will log to the console
## logfiles can be set with 'set_logfiles' and 'add_logfiles'

set_logfiles(logger, #logger can be an object or a character string
             c("test.log"))
logger$m_files
jlog.info(logger, "Test") ## Output to "test.log"
add_logfiles(logger, #DITTO
             "")
logger$m_files
jlog.info(logger, "Test 2") ## Output to both "test.log" and the console

## The logger will automatically append the new messages to the files.
## logfiles can be flushed with 'jlflush' function
jlflush(logger) #Flushing console isn't relevant, 'test.log' will however be flushed

## Prefix
## The prefix is the string printed right after the date. It can be set directly through the member logger$m_prefix
logger$m_prefix = "Michael Jackson"
jlog.info(logger, "is bad")

## Logger can be fully configured at inception. 'JLoggerFactory ' accepts additional arguments.
## Namely 'files', 'prefix' defaulting to the logger's name, 'level'
my_logger <- JLoggerFactory(name = "my logger",
                            files = c("", "file1.log", "file2.log"),
                            prefix = "make it funky",
                            level = JLOGGER.TRACE)
my_logger
jlog.trace(my_logger, "yaay!")

## The default level is JLOGGER.DEFAULT_LEVEL
JLOGGER.DEFAULT_LEVEL
## If changed the next jlogger instantiated without a specified level will be set to this level
JLOGGER.DEFAULT_LEVEL = JLOGGER.ERROR
errlog = JLoggerFactory(name = "err logger")
errlog$m_level

## Turning off logging completely
## Logging can be turned off by setting a logging level higher than the highest logging level available. The highest level for now is JLOGGER.FATAL
logger$m_level = JLOGGER.FATAL + 1L
jlog.fatal(logger, "1 = 2") # won't print

## The jlog functions accept NULL as input. In this case no message is printed

jlog.info(NULL, "1 != 1") ## Nothing
## It makes it easy to turn off all messaging without having to handle logging levels.


## Clearing registered loggers.
## You can clear all created logger with 'JLoggerReset'
JLoggerReset() ## subsequent calls to JLoggerFactory will create a new instance

## jlogger config
## JLoggers can be configured through a config list
## Logger configs are used to automatically configure the logging leveld
## You need to specify one entry by logging level. You don't need to specify all the levels. Each entry refers to a set of regular expressions.
## Example:

lconf <- list(ERROR = "err.*",
              INFO = c("jalgos", "inf.*"),
              TRACE = c("tracer", "overlogger"))

## Now by setting the 'logconfig' argument in JLoggerFactory to 'lconf', the name will be matched against the different regular expression and in case there is a match the logging level will be set to the corresponding level.
## If a name doesn't match any regex in the config, the logger will be quiet for all levels

err_log <- JLoggerFactory(name = "err logger",
                          logconfig = lconf)

err_log ## "err logger"  matches "err.*" hence err_log's level is set to ERROR

jalgos_logger <- JLoggerFactory(name = "jalgos",
                          logconfig = lconf)
jalgos_logger ## "jalgos" matches an entry for INFO
