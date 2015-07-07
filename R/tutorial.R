## This tutorial aims at introducing the R language to programmers who are already familiar with other languages.

## More resource is available here: http://cran.r-project.org/
## For a little history: https://en.wikipedia.org/wiki/R_(programming_language)

## The main attractive aspect of R is that it is a high-level language whose syntax remains close from mathematical approximation and leaves out the technical aspects of implementation

## R has an invite (`R` in linux) where the user can enter a sequence of commands and see the results in the invite. This speeds up development considerably as results can easily be visualized on the fly.
## R is interpreted and does not require compulation, it can be use as a script with the binary `RScript`

## R language is functional meaning that a call f(x) will not modify x, though there a few exceptions.
## The language is very flexible and is nowhere far from the rigor of other languages like C++ or Java

## R was designed for statistical analysis and relies on wide collection of package than can perform low level mathematical computation.
## Since R is open source, it has a huge number of users and contributors who constantly submit new packages to perform analysis on very various domain.

## R is used by writing the higher logic of your analysis in R and calling packages functions to perform the costly computations.

## R is sometimes deemed as a "slow" language but a well written R program should have similar performance to its C/ C++ counterpart and development should be way faster in R

################## Tutorial #######################

## Let's start
##################################################
## I) Basic functionalities                      
##################################################


## First open an R invite by typing R in a terminal
print("Hello world!") # print is the function to print

## Basic operation:
1 + 1 # should print 2

## There a few basic R types. notably
1L ## integer range [-2 ^ 31 +1; 2 ^ 31 -1]
1.5 ## numeric
TRUE ## logical
"a" ## character

## R also handles complex and other types
1 + 1i ## is complex

## It is important to note that everything is a vector in R
## R was design with the idea to solve scientific problems hence vectorization made sense
is.vector(1) ## TRUE

## Assigning a variable is trivial:
x = 1
## R has two operators for assignatio `=` and `<-`. So above is equivalent to:
x <- 1

## R syntax doesn't require any final semicolon an end of lineis interpreted as the end of a command. Semicolons can separate a sequence of command on the same line
a = 1; a = a + 1; a * 2

## R language offers the classical imperative features found in other languages: if else, for/ while loops
if(1 == 2){
    print("No!")
} else {
    print("Yay!")
}

while(x != 10){
    print(x)
    x = x + 1
}

## iterating through a vector is done like this
my_vector = 1:10 #`:` is an operator to create a sequence from n1 to n2
for(x in my_vector){
    print(x)
}

## !!!! Important note on loops !!!!!
## A good rule of thumbs in R is "avoid using loops". Indeed since R is interpreted an R loop is painfully slow
## The real thing to avoid is looping through your data
## For loops can be used if what is performed during the loop is what takes the most time
## R spirit's is to externalize all the big computatio and only write the higher logic
## e.g: sum of a vector
r = rnorm(1000000) #huge vector of n random draws of a N(0, 1)
system.time({x=0; for(i in 1:length(r)) x = x + r[i]; print(x)}) ##slow
system.time(print(sum(x))) ## lightning fast

##################################################
## II) Functions
##################################################

## Defining functions is straightforward
my_function <- function(x, y)
{
    x + y - y / x
}

## calling my_function:

my_function(exp(1), pi) ## pi is Pi 3.14...

## functions can be called with named arguments

my_function(y = exp(1), x = pi) ##order don't matter names matter

## functions can have default values

other_function <- function(x,
                           p = 10)
{
    x ^ p
}

other_function(2)

## ellipsis
## R has a way to pass down extra arguments to other functions, the ubiquitous and infamous ellipsis: `...`

my_function_with_ellipsis <- function(x,
                                      ...)
{
    other_function(x = my_function(x, ...), ...) #Note that ther is no conflict between the name x and the variable
}

my_function_with_ellipsis(x, 2) #equivalent to  other_function(x = my_function(x, 2), 2)

## it is really useful to pass down configuration parameters.
