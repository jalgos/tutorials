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

## Help

## To get help on any function use `?`
?cat ## opens the manual of the cat function
## To look for a word in the entire R manual use `??`
??cholesky ## Will lookup for all the occurences of the word cholesky in the help


##################################################
## III) Basic data structures
##################################################
## Vectors


##################################################
## III) Functions
##################################################

## Defining functions in R is pretty simple. Function definition offers a lot of flexibility to implement the high-level logic of the data analysis.
## Types of arguments are not set beforehand. The execution of a fonction will end if the type of an argument is not handle by an operation which makes R functions inherently generic. 
## 

## Defining functions is straightforward
my_function <- function(x, y)
{
    x + y - y / x
}

## The above function will work on numeric, logical, complex vectors or any data type that has + - / defined

## Functions will not modify its arguments (passing by value instead of passing by reference)
## Functional languages have many advantages. The two main advantages are that termination of an algorithm can be proved and debugging is made infinitely easier as objects are stateless.
mod_func <- function(x)
{
    x = 2 ^ x
    x * (1 - x)
}

x = 3
mod_func(x)
x ## still 3

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
## `...` is a sequence of extra arguments. `...` can contain any variable type and can be any arbitrary size.

my_function_with_ellipsis <- function(x,
                                      ...)
{
    other_function(x = my_function(x, ...), ...) #Note that ther is no conflict between the name x and the variable
}

my_function_with_ellipsis(x, 2) #equivalent to  other_function(x = my_function(x, 2), 2)

## it is really useful to pass down configuration parameters.
## If arguments are not named they will be taken in the order they appear.
a = 2
b = 9
other_function(a, b) ## equivalent to other_function(x = a, p = b)
## Naming the arguments makes up for a less ambiguous code. And it greatly helps passing down arguments to other functions
## If a named argument is not in the function's signature, R will yield an error. If you give more arguments than in the function's definition R will stop with an error too

cos(p = 10) ## error the function is defined for an argument named 'x'
cos(2, 3)  ## error cos takes only one argument

## This happens unless the function takes the ellipsis in its definition

fune <- function(x, ...) x + 8

fune(1, 10, list(20), "A", h = "OOOO") ## no error

## It's also possible to define a function that takes an arbitrary number of extra arguments. This is similar to C++ extra arguments feature besides that R does not put any requirement on the type of the extra arguments.

## Iterating through `...` is done with creating a list containing the element of `...` and iterating through them

print_class <- function(...)
{
    for(item in list(...))
    {
        print(class(item))
    }
}

print_class(1, 1L, TRUE, list(a = 1), NULL, "A")

## There is another way to access the nth element of `...`: `..n`
print_2class <- function(...)
{
    print(class(..1))
    print(class(..2))
}

## However this feature is only here for debug purposes and should not be used in your function definition, named arguments are to be used instead.

## The position of the `...` in the function matters a lot.
fab <- function(a = 1, b = 2, ...) a + b
fab(2, 3) ## 5

fxy <- function(x = 1, ..., y = 2) x + y
fxy(2, 3) ##4

## In the latter case, the second argument was interpreted as part of `...`. The only way to override y is by explicitly naming it
fxy(2, y = 3) ##5

## Global assignment.
## It is possible to globally assign a function from a function thanks to the operator: `<<-`

fa <- function(y)
{
    MYGLOBALVAR <<- y
}

fa(2)
MYGLOBALVAR

## Global assignment are to be avoided. It is usually a bad idea to try to circumvent R's functional behavior. Configuration variables can be passed down to functions.

## Functions can be arguments of other functions

fap <- function(x, f) f(x)
fap(3, exp)

## The missing function
## To test if an argument was passed to a function, one can use the function missing

fmiss <- function(x, y)
{
    if(missing(x)) print("x is missing")
    if(missing(y)) print("y is missing")    
}

fmiss(1, 2)
fmiss(1)
fmiss(y = 1)




##################################################
## IV) R Objects and method dispatching
##################################################

## Method dispatching

## We saw earlier that regular R function can already be very generic, but there is a proper way to achieve polymorphism thanks to generic function

## Let say we want to define a function plus that concatenate two strings for variables of type character and performs an addition for variables of type numeric

## First we define a generic
setGeneric(name = "plus", def = function(x, y) standardGeneric("plus"))

## The call above a generic (polymorphic) function that takes two arguments named x and y. Technically it adds an entry in the mtable, similarly to what happens when we define a C++ object.
## The use of standardGeneric means that there are now concrete default making "plus" virtual. It is possible to give a default function.
## e.g:
## setGeneric("plus", def = function(x, y) x + y))
    
plus.number <- function(x, y) x + y
plus.character <- function(x, y) paste0(x, y)

## Now concrete implementations for set types can be added with the setMethod function
## the arguments are the generic function name, the types of the arguments to dispatch the function on, and the function definition.

setMethod(f = "plus", signature = c(x = "numeric", y = "numeric"), definition = plus.number)
setMethod(f = "plus", signature = c(x = "character", y = "character"), definition = plus.character)
plus("A", "B") ## "AB"
plus(2, 1) ## 3

## R provides a generic data type: ANY. It matches any data type
## If we want plus("A", 1) to return "A1" we can do this:

setMethod(f = "plus", signature = c(x = "character", y = "ANY"), definition = function(x, y) plus(x, as.character(y)))
setMethod(f = "plus", signature = c(x = "ANY", y = "character"), definition = function(x, y) plus(as.character(x), y)) ## The symmetric
plus(1, "A")
plus(1, 2)
plus("A", 1)

## The key word 'missing' can be used in the signature.
setMethod(f = "plus", signature = c(x = "ANY", y = "missing"), function(x, y) plus("bob", x))
setMethod(f = "plus", signature = c(x = "missing", y = "ANY"), function(x, y) plus(y, 10))
plus(1)
plus(y = 1)


## Defining operators
setGeneric(name = "%+%", function(x, y) plus(x, y))
setMethod("%+%", c("ANY", "ANY"), plus)



## R Objects

## This is a more advanced feature of R
## There are 3 types of R objects: S3, S4 and R5

## S4 is the successor of S3 and are a way of defining functional objects or traits.
## R5 are mutable objects based on S4. We strongly advise you to avoid using them.

## We will study the S4 syntax and we recommand that you use S4 over S3 anytime you define a new object.

## S4 objects are functional objects. It means that if X is an S4 object X won't be modified by the call f(X). 
## Unlike Java or C++ objects, R objects are immutable.
## S4 object don't have dedicated method. Polymorphism is achieved through R's function dispatching.
## S4 objects are used to achieve encapsulation and inheritance.
## As in other languages S4 objects can be abstract and can inherent from one another. Multiple inheritence is allowd.
## S4 classes are defined with the setClass function
A <- setClass("A",
              slots = c(x1 = "numeric", x2 = "numeric", name = "character"))

## new instances can be created by a call to new:
myA = new("A", x1 = 10, x2 = 5, name = "atest")

## setClass returns a function that is equivalent to new("A", ...)
## So above call is equivalent to:
myA = A(x1 = 10, x2 = 5, name = "atest")

##Object components are accessed with `@`
myA@x1

## Feel free to redefine the constructor as more complex objects require a finer way to be initialized

## inheritance is achieved by the 'contains' argument in setClass

B <- setClass("B",
              slots = c(x3 = "numeric"),
              contains = "A") ## multiple inheritance is allowed
myB1 = new("B", x1 = -9, x2 = -5, x3 = pi, name = "btest1")

## The constructor new will assume all named parameters are class members and unnamed parameters are superobjects.
## Following object is build on top of myA
myB2 = new("B", x3 = 2.3, myA)

## Inheritance from the class "VIRTUAL" makes the class abstract and not instantiable

V <- setClass("V",
              contains = "VIRTUAL",
              slots = c(u = "logical"))

V(u = TRUE) ## Error cannot instantiate

C <- setClass("C",
              contains = c("V", "A"), #multiple inheritance
              slots = c(z = "ANY"))

myC = new("C", u = FALSE, z = list(a = 1), x1 = 0, x2 = 8, name = "ctest") ## Can be instantiated

## The function is can determine whether an object is of a given class

is(myA, "A") ## of course
is(myB, "A") ## inherits
is(myA, "B") ## B is a subclass of A so no
is(myB, "B") ## obviously
is(myC, "V") ## inherits
is(myC, "B") ## same super class A but different

## Writing functions
myfA <- function()
