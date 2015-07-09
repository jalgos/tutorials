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
## R spirit's is to externalize all the big computation and only write the higher logic
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
## II) Basic data structures
##################################################
## Vectors

## In R all basic types are vector
is.vector(1) ## TRUE
is.vector("A") ## TRUE too
## As a statistical language R was designed to analyze collections, sequences, tables of data. So vectorizing the basic types makes sense.
## This design feature makes vector computation faster.

## It is really important to keep this design in mind when programming in R

## Function must be design to work with vectors and not individual elements.

## All the basic operators apply term by term
## e.g if we have two vectors X = (x1, ..., xn) and Y = (y1, ..., yn)
## Then X ^ Y = (x1 ^ y1, ..., xn ^ yn)
X = rnorm(100) ## A vector of length 100 with each element drawn from a normal 0, 1
Y = rnorm(100) ## DITTO
X + Y
X > 0 ## A vector of type 'logical'
cos(X)

##Creating a vector
## The function 'c' is use to concatenate vectors
A = c(1, 2, 4)

## 'length' returns the size of a vector
length(A)
length(X)

## Subsetting vectors
## A subset of a vector can be extracted with the`[]` operator
## Its arguments are a vector of indices (numeric or double) or a vector of booleans equal to the length of the vector to be subset
X[A] ## keeping indices 1, 2, 4

## !!! Important !!! R indexing starts at 1 not 0
X[0]

X[Y > 0] ## Subset of X for which Y is positive
X[] ## no parameters provided, returns the whole vector
## indices can be repeated. If an index exceeds the size an NA will be returned

X[c(2, 2, 1, 1000)]

## Negative indices can be use to return the complement of the subset formed by the indices given
X[-(11:80)] ##Will return the first 10 cells and the last 20

## Creating a vector from 1 to n is easy:
1:10

## The function 'seq' can be used too
seq(from = 1, to = 10, by = 1)
## seq is not limited to integer
seq(from = 0, to = 1, by = .1)

## There are a lot of native function to aggregate vectors
sum(X)
max(X)
mean(X)
quantile(X, c(0.1, 0.2, 0.8))

## names

## Vector can have name attributes
## To assign names to a vector simply use the function 'names<-'
u = rnorm(26)
names(u) = letters
u
## names must be characters
## If u is named u can be subset by names
u[c("a", "l", "m")]

## Names can be set when explicitly creating a vector
V = c(u = 1, bob = 2, x = 9, 10, 22, f = 30)
names(V) ## not every cell needs to be named


## Matrices

## R has a native type to perform linear algebra

M = matrix(rnorm(100), 10, 10) ##10 x 10 matrix
## dimensions can be queried with dim, nrow, ncol
dim(M)
ncol(M)
nrow(M)
## '%*%' is the matrix multiplication operator
M %*% t(M)

## Multiplying a matrix by a vector returns a vector
V = rnorm(10)
M %*% V

## Matrices can be subset similarly to vectors
M[1:4, 5:8]
## if a subset index is missing no subsetting is done on that dimension
M[1:4, ]
M[, 5:8]
M[, ] ## Whole matrix

## R will convert the result to a vector if any of the subsetting dimension is of length 1
M[5, ]
M[, 10]

## This can be avoided by setting drop to FALSE
M[5, , drop = FALSE]
M[, 10, drop = FALSE]

## Negative subsetting works with matrices too
M[-(1:5), 1:5]

## Matrices can be converted to vectors
V = as.vector(M)
## M[i, j] = V[i + (j - 1) * n] where n is the number of rows in M
## This is also known as column major representation
## Hence M can be subset as a vector:

M[18] == M[8, 2] ## TRUE

## Similarly to vectors matrices can have row names and column names
## The functions use to access and set names are 'rownames' 'colnames' and 'dimnames'
rownames(M) = letters[1:10]
colnames(M) = letters[11:20]
dimnames(M) ## list of rownames and colnames
M[c("a", "b"), c("m", "p")]
## Names are elegantly assigned when performing most of the operations on matrices
M %*% t(M) ## rownames and colnames are identical now

## R has plenty of native functions to analyse matrices:
det(M)
solve(M) ## inverse
eigen(M) ## eigen value decomposition

## R supports sparse matrices thanks to the package 'Matrix'.
## It is installed by default.
## Sparse Matrices are used in many kinds of problems, notably in graph analysis.





## List
## R provides a very flexible list structure
## R list can be composed of any kind of objects. 
L = list(a = 1, b = "bonjour", r = rnorm(1000), l = list(1, 2, 3, 4), 12, cos, fun = exp)

## Similarly to vectors list have names
names(L)
## list are access with the `[[]]` operator
L[[1]]
L[["a"]]
## The `$` operator is another way to access members
L$l

## The `[]` returns a list
L[c(1, 2, 4)]
L[c("a", "fun")]

## Lists are a great way to compose more complex objects without having to resort to writing a class.



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

## By default the return value is the last evaluated call
## R has a return statement that can be use to return early

my_function <- function(x, y)
{
    if(any(!is.finite(x))) return(0)
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
## Binary operators can be define/ overriden. There is a finite set of binary operators allowed in R. '%+%' is one of them
setGeneric(name = "%+%", function(x, y) plus(x, y))
1 %+% "A"
2 %+% 1

## This feature is useful for overriding mathematical operators for customly defined objects (e.g. matrix multiplicatiob...)

## R provides a quick way to define a virtual superclass of a define set of classes. This is a useful feature to define a function that can handle several different types.

setClassUnion("allnumeric", c("numeric", "complex", "logical", "integer"))
setMethod(f = "plus", signature = c(x = "allnumeric", y = "allnumeric"), definition = plus.number)
1L %+% 2 + 0.5i

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
myfA.A <- function(X)
{
    print(X@name)
}

setGeneric("myfA", function(X) standardGeneric("myfA"))
setMethod("myfA", "A", myfA.A)

myfA.C <- function(X)
{
    print(X@u)
}

setMethod("myfA", "C", myfA.C)
myfA(myA)
myfA(myB1)
myfA(myC)

## A great advantage of function dispatching is that it can be done on several arguments.
setGeneric("combineA", function(X, Y) standardGeneric("combineA"))
setMethod("combineA", c("A", "A"), function(X, Y) "AA")
setMethod("combineA", c("B", "B"), function(X, Y) "BB")
setMethod("combineA", c("A", "C"), function(X, Y) "AC")
setMethod("combineA", c("B", "C"), function(X, Y) "BC")

combineA(myA, myA) ## AA
combineA(myB1, myB2) ## BB
combineA(myB1, myC) ## BC
combineA(myC, myB1) ## AA
combineA(myC, myC) ## AC

## Dispatching functions on multiple arguments allows for finer control

## "Modifying" objects

## As mentionned earlier R objects don't hold states. Modifiying an object is done functionnaly

modify <- function(X, new_name)
{
    X@name = new_name
    X
}

newA = modify(myA, "newA")
## It is preferable to use different variables instead of reassigning the value of modify to myA.
## Both syntax are correct but assigning the output of modify to a new objects makes it easier to debug.
## Of course if myA is a huge object, then it's better not to systematically copy it.

## Defining objects is an advance feature that is usually used at an advanced stage in the research. Basic exploratory analysis do not require objects as structures like lists, data frames and vectors can easily do the job.
## It is however important to know how to implement R objects to better understand some R packages and to be able to soundly design algorithms.


##################################################
## V) Debugging functions
##################################################

## R has a step by step debugger somehow a bit more basic than the C++ or Java equivalents but yet very powerful

## A breakpoint on the entry of a function can be set with the debug function
fun_test <- function(x)
{
    x = x + 1
    for(i in 1:10)
    {
        x = x * 1.1
    }
    x
}
debug(fun_test)
fun_test(1)

## You've entered the debugger. To go to the next line enter 'n'
## To continue execution until the next breakpoint type 'c' or 'cont'.
## To see the call stack enter 'where'
## 'Q' to quit the browser
## Typing 'c' from within a loop will continue execution until the end of the loop.
## Try the different debugging commands with fun_test

## To stop debugging a function use 'undebug'
undebug(fun_test)
fun_test(1)

## If you only want to debug a function and not debug it at the next call use 'debugonce'
debugonce(fun_test)
fun_test(1)
fun_test(1)

## Another way to invoke the debugger is calling the function browser
browser() ## No real use from the top level

## If inserted in a function's body it will stop there
fun_test <- function(x)
{
    x = x + 1
    for(i in 1:10)
    {
        x = x * 1.1
    }
    browser()
    x
}
fun_test(1)

## There is another way to set breakpoints at a specific line of a file or of a function: 'setBreakpoint' and 'trace(..., at = x)'
## Most of the time debug will cover most of your needs.
## Please refer to the manual if you want to use these functions

## Debugging generics

## If functions are dispatched on their arguments, 'debug' will not do what you want it to do.
## You need to use 'trace' and specify
setGeneric("ftest", function(x) standardGeneric("ftest"))
setMethod("ftest", c(x = "numeric"), function(x)
          {
              print("BLAH")
              y = 10 * x
              cos(y)
          })
setMethod("ftest", c(x = "integer"), function(x)
          {
              print("BLUH")
              y = x / 10
              exp(y)
          })
debugonce(ftest)
ftest(1)
## trace is a finer version of debug
## You can set the function it invokes when hitting the breakpoint, specify behavior when exitting the function, stop at a specific line within the function.
## The feature we'll look at is selecting the function to trace according to its signature

trace(ftest, signature = "numeric", tracer = browser)
ftest(10) ## stops
ftest(1L) ## integer version is not debugged

## use untrace to remove breakpoints
untrace(ftest, signature = "numeric")

## Debugging is fondamental in developpment. It is really important to master this tool to develop faster in R

##################################################
## VI) Using packages
##################################################

## R provides a myriad of packages written in faster languages such as C/ C++ or Fortran.
## It is strongly advised to use packages instead of writing your own structures.
## Installing a package is done with the 'install.packages' command
install.packages('data.table') ## Great package that extends R's dat frames
## To load a package use the 'library' command
library(data.table)

## This is pretty straightforward
## Before starting to implement something that is not central in your project make sure to thoroughly look on CRAN if a package is available for what you want to do.
