## The following tutorial, written by Sebastien Lamy, data scientist at Jalgos, is meant to be executed line by line.
## It will give a good overview of the data table library for R, a highly efficient tool to do quick and clean data exploration and analysis.
## This tutorial supposes a basic knowledge of R and dataframes.
## The comments will walk through each step.
## Enjoy!

## uses the ESPNFC datasets on football

############## DATA TABLE INTRODUCTION ###################

## Let's play around with the ESPNFC data set to understand how datatable works

## To load the data and sources and see a few example visualizations, run the following command in an R shell opened at the root of the project (../R_tuto/) and press keys until you are back on a terminal prompt.

source("analysis/example_analysis.R")

## A data table extends a dataframe so it has all the dataframe's functionality and properties
## In particular it has a [] operator to query its content
## Its implementation is simple. It's a list of same lengths vectors. The vectors are the columns.
## To get help on data.table type:
## ?data.table
## the [] operator has several arguments. The most important are the first two arguments i and j.
## i is used to apply a filter on the rows
## j is used to work on the columns
## We'll start working on a data table object named GOALS.

## All the arguments are optional so the following call will return the entire set

GOALS[]

##### Use of i
## The most simple use of i is to give it a vector of indices to get a subset of the data

## this will return the first line of the dataset (indexing starts at 1 in R, not 0)

GOALS[1]

## this will return a sub dataset with the 1st line, 3rd line twice and the 5000th line

GOALS[c(1, 3, 3, 5000)]

## An alternative is to give a boolean vector which size is the same as the table's number of rows
## How to get the number of rows of data.table/frame

nrows = dim(GOALS)[1]
## dim(GOALS) returns a vector of 2 values: number of rows and number of column. Once again, be careful, vectors indexing start at 1 in R!!

## The following instruction takes a random sample of GOALS

GOALS[sample(c(TRUE, FALSE), nrows, replace = TRUE)]

## A little explanation here: sample(vec,size,replace=TRUE) returns a vector of length "size" containing a random sample of the values contained in "vec". So we just made a simple call working on data table "i" parameter to filter on the lines.

## A great feature of data.table is that the columns are easily accessible from within the data.table
## First, let's have a look at our columns:
colnames(GOALS)

## We see that GOALS has a column goal_number. It can be directly accessed in j this way:

## i is left missing, no filtering will be done on the lines.
GOALS[, goal_number]

## this is equivalent to:
GOALS$goal_number
GOALS[["goal_number"]]
GOALS[j = goal_number]

## So now if we want to filter GOALS on the lines where goal_number equal to 6:
GOALS[goal_number == 6]

## If we want only penalties:
GOALS[penal == "pen"]

## i can be any function that evaluates to an integer vector or to a boolean vector of the right size (th right size being the number of lines of your data table)
## e.g all goals where the square of the goal_number is greater than the minute
GOALS[goal_number ^ 2 > minute]

## less silly (or is it?): all goals scored by a van something.
## If you are familiar with regex this should be crystal clear:
GOALS[grep("\\svan\\s", scorer, ignore.case = TRUE)]

###### use of j
## We saw earlier that we can get the content of a column as a vector like this:
GOALS[, scorer]

## You can get a subset of the columns like this:
GOALS[, list(scorer, player_id, team_id)]
## The result is another datatable!

## if we want to get a subset and rename the columns
GOALS[,list(sc = scorer, pid = player_id, tid = team_id)]


## j will evaluate expression looking up variables in the datatable then in the toplevel environment
x = c("AHA", "IHI")
## no columns AHA or IHI in GOALS, so check out what this call returns:
GOALS[, x]
## AHA IHI

## It's really important to remember that in j columns become named variables. So you can't access them through strings or index.

GOALS[, "scorer"]## scorer string not the scorer column
GOALS[, 1]## constant 1, not the first column

## If you want to access the columns of a datatable without knowing beforehand the column names (it happens a lot when doing advanced stuff) there are several ways that we will comment later.

## The good thing about j is that you can enter any function  in it, and it will just be executed. You can for example plot things directly inside j.

GOALS[, plot(goal_number,minute)]

## sum/average columns
GOALS[, mean(minute)]

####### Adding/ removing columns
## to add a column or remove a column is simple with the := operator
## adding
GOALS[, super_column := goal_number + minute]## this won't output any result
GOALS

## removing is even easier
GOALS[, super_column := NULL]
GOALS

## To assign multiple columns at the same time:
GOALS[, c("col1","col2") := list(goal_number / 2, cos(minute))]
GOALS

## To remove several columns at the same time:
GOALS[, c("col1","col2") := NULL]
GOALS


######### the "by" argument
## A key feature in data.table is the ability to apply function by groups

## e.g median minute by goal type
GOALS[, median(minute), by = penal]

## the result is a datatable itself. The unnamed columns will be V1...Vn
##  to name them just do
GOALS[,list(med = median(minute), max = max(minute)),by = penal]

## There are a few hidden variables in the datatable.
## One of them is .N It gives the number of rows
GOALS[, .N]
## The number of goals by player can be obtain as such
NG = GOALS[, list(num_goals = .N), by = scorer]
## ordered and keep the top 25
NG[order(num_goals, decreasing = TRUE)][1:25]



######### A short explanation about the syntax
## The order of arguments in [] is i,j,by,...(other arguments that we will explain later). To use a similar logic to dataframe's we omitted the names of the argument so far.
## That's why we like not to name i and j.
## One could also ommit naming the by argument. The above expression is equivalent to:
NG = GOALS[, list(num_goals =.N), scorer]## Since by is the 3rd argument of []
## One could also order the arguments the way he wants by naming them
NG = GOALS[by = scorer, i = minute > 45, j = list(num_goals = .N)]## number of goals by player in  the second half

## But it's better to follow the convention and not name i and j and name the rest of the arguments. Naming i and j is acceptable if it is confusing to you. But not naming any arguments is to be avoided.



####### Keys
## A very important feature of data.table is the ability to index a dataset by a set of keys to be able to make faster queries. Keys are also critical in joining tables
## The keys combination doesn't have to be unique
## example:
setkey(GOALS, team_id)
## equivalent to:
setkeyv(GOALS, "team_id")## pay attention to the "v"
## The order of GOALS has changed:
GOALS
GOALS["125"]## This will get all the goals scored by team_id 125

## The data can be fully indexed by game_id,team_id,goal_number
setkey(GOALS, game_id, team_id, goal_number)
GOALS[list("99019", "448", 2)]
## check the keys of a data.table
key(GOALS)
## shortcut to set all columns as keys
setkey(GOALS)
## unset keys
setkey(GOALS, NULL)

## the "keyby" argument. "keyby" can be used instead of the "by" argument to return a keyed datatable by the groups used.

NG = GOALS[, list(num_goals = .N), keyby = list(scorer, team_id)]
key(NG)


##### Joining tables

## Joins are essential when manipulating data. In this particular project we have several tables containing informations about teams, leagues, players. If we want to add the information of a table into another one we do it through joins.
## People who are familiar with SQL will understand the concept and find it way easier to grasp and more minimal
## As an example we'll get the number of goals by league
## The league column is not in the goal table. We'll add it by joining GOALS with the games table GS
## The first way to join to tables is through i
## We can give a datatable in the "i" argument as this: D1[D2]. D1 and D2 must be keyed and the intersection of their keyset must be non empty. The result will be a new data.table composed of the union of the records of D2 with their matches in D1. It is a "one to many" join in the sense that if there are several matches in D1 for a record in D2 all the combinations will appear.
## Let's work on a trivial example

D1 = data.table(a = c("A", "A", "B", "B", "D", "E"), b = (1:6) ^ 2, name = c("JIM", "JON", "ART", "MO", "BO", "C.J"))

D2 = data.table(a = c("A", "O", "Z", "B", "E"), y = rnorm(5),
    name = c("JIM", "KIM", "DAN", "ART", "A.J"))

## let's do a join on the names. We need to key both tables on names
setkey(D1, name)
setkey(D2, name)

## Let's do the join
D1[D2]## As many records as in D2. The records that don't have any match in D1 are filled with NA values. Duplicate columns will have a .1 suffix added to them. The column with the .1 is the column of the data.table in the "i" argument

D2[D1]## Pay good attention to the difference!

## Both the above joins where 1 to 1. let's do a one to many join. The a column has repeated values.

setkey(D1, a)
setkey(D2, a)
D2[D1]## still 1 to 1
D1[D2]## one to many. JIM and ART in D2 are matched to several records in D1
## The behavior of a one to many join can be controlled by the parameters allow.cartesian and mult
## if allow.cartesian is set to FALSE, a one to many join will fail. This is useful if you know beforehand that the join should be one to one and you want to catch errors
D1[D2, allow.cartesian = FALSE]
## The default behaviour is controlled by the option "datatable.allow.cartesian" c.f ?data.table. Here is how it's set:
## options(datatable.allow.cartesian=TRUE)
## mult will tell which records to match in case of a one to many join. It can be either "all" (default), "first" or "last"

D1[D2, mult = "first"]
D1[D2, mult = "last"]


## Works also with several keys

setkey(D1, name, a)## name and a have to match
setkey(D2, name, a)
D1[D2]
D2[D1]

## practical example
## let's join GS and GOALS on the game_id column
setkey(GOALS, game_id)
setkey(GS, game_id)
GS[GOALS]## A huge table! It's better not to assign it to anything to save memory and work directly on it
GBL = GS[GOALS][,list(num_goals = .N),by = league_id]
## We want to add the league name now. It's in LEAGUES
setkey(GBL, league_id)## note we could have used keyby to create the table and save us this call
setkey(LEAGUES, league_id)
GBL = LEAGUES[GBL]
GBL[order(num_goals, decreasing = TRUE)][1:25]

####### using joins to add columns.
## We saw the use of operator ":=" earlier to add a new column. This can be done when joining too
## e.g instead of producing a huge table by joining GS and GOALS we can add the league_id column to GOALS. This is done this way:

GOALS[GS, league_id := league_id]## This can seem confusing as league_id is on both side but data.table understand the LHS of := is a new column and will look up the LHS using the usual rules

####### merges
## So far joins were done in an asymetrical way. D1 and D2 play a different role. D1[D2] will just add D1's info to D2

## The "merge" function provides a way to join tables in a way where D1 and D2 play symmetrical roles
## In the following example will merge to kind of player stats. One obtained from GOALS and the other from BXS

NG = GOALS[, list(num_goals = .N), by = list(player_id)]
NAO = BXS[is.finite(a) | is.finite(of), list(assists = sum(a, na.rm = TRUE), offsides = sum(of, na.rm = TRUE)), by = player_id]

## Both tables contain information of equal importance. So we want to keep all the info in both table when merging. This is done like this

merge(NG, NAO, all =TRUE, by = "player_id")## all means keep all records in both table. The by argument is a string and it doesn't matter how the tables are keyed beforehand.
## You can do an "asymmetrical join" by specifying all.x or all.y. Setting all will set all.x and all.y to the same value

merge(NG, NAO, all.x = TRUE, by = "player_id")## Discard non matching values in NAO. Equivalent to NAO[NG] if keys were set properly.

merge(NG, NAO, all.y = TRUE, by = "player_id")## equivalent to merge(NAO,NG,all.x = TRUE,by = "player_id")

## If we want to keep the intersection of the matching keys then we set all = FALSE (or don't set anything at all as it is the default behaviour)

merge(NG, NAO, by = "player_id")


#######rolls
## Rolls are used to match records to the nearest value. This is useful for example when dealing with time series. If we want to get the latest value available before a date we will use roll.

## Example. Assume we make a simple model were we consider that players are a certain "state", characterized by one parameter, that can change every 5 minutes during a game. We will first create the table that will contain this parameter: we need to have, for each player and each game, one line every five minutes of the game. So we will make a cartesian product of the list of players (through player_id) and of the 5 minutes sets that a 90 minutes game contains. We thus have a table of the good size to contain the parameter we want:

PARS = cartesian_data_table(NG[, list(player_id)], data.table(minute = seq(0, 90, 5)))

## Let's fill it up with a parameter. We start with a randomly generated parameter, using function rnorm (random generation of the gaussian function) with the argument .N in j to have the proper number of values.
PARS[, par := rnorm(.N)]
## Let's have a look:
PARS

## PARS is a table of parameters. Each player has a parameter for every five minute.
## To fit our model we want to get the parameter that the player has when he scores. Goals are scored every minute, so we need to get the parameter for the closest 5 minute bucket in the PARS table.
## This can be done through a roll. We'll joing GOALS and PARS through player_id and minute. Let's set the respective keys:

setkey(PARS, player_id, minute)
setkey(GOALS, player_id, minute)
## The roll can only be done on the last key, that's why it's important that minute is the last key.

PARS[GOALS]## is not satisfactory a lot of record won't have a parameter

PARS[GOALS, roll = TRUE]## With roll=TRUE parameters will be match to the GOALS where minute is just greater than minute in PARS
PARS[GOALS, roll = "nearest"]## This will match to the nearest minute

## Be aware that PARS[GOALS,roll=TRUE] and GOALS[PARS,roll=TRUE] will have a different result. It can ruin an analysis if precedence is critical (e.g. finance, prediction)

## We'll add the parameters column this way

GOALS[PARS, par := par, roll = "nearest"]


##### Dealing with data.table programmatically
## In many cases we don't want to explicitly write the column names in the expressions. Because either we don't know them or because we want to apply the same function to all the columns.
## Remember that a datatable is a dataframe and hence a list
## You can extract columns by their names (character vector) in multiple ways
GOALS[["scorer"]]## not really useful as you lose all the data.table features
GOALS[c("scorer", "player_id")]## Won't do what you want it to do
as.data.frame(GOALS)[c("scorer", "player_id")]## Works but ugly and slow

## data.table provides a way to extract columns by their names with the argument "with". If set to false, it will assume that j is a list of column names
GOALS[, c("scorer","player_id"), with = FALSE]

## There is another hidden variable: .SD. It contains the list of the columns.
GOALS[, .SD]## simply itself

## Note that if we assign the result of this call to a variable, the resulting data.table will be read only

TEST = GOALS[, .SD]
TEST[, a := 1]
TEST[, scorer := NULL]

## You can use the .SDcols argument to limit the columns you want to study

GOALS[, .SD, .SDcols = c(1, 2, 3, 4, 5)]
## use match to get the columns positions (check ?match if you don't know it)
col_pos = match(c("scorer","team_id"), names(GOALS))
GOALS[, .SD, .SDcols = col_pos]

## .SD is useful when we want to apply a function to a subset of the columns

GOALS[, lapply(.SD, mean, na.rm = TRUE), .SDcols = match(c("goal_number", "extra_time", "minute", "par"), names(GOALS)), by = penal]

## Another way to handle data programmatically (without knowing the columns beforehand) is to use expressions

## expressions are a R structure. It is not evaluated unless eval is called upon them
exexp = expression(5+1)
eval(exexp)
x = 10
exexp2 = expression(x^2)
eval(exexp2)
## let's say we want to eval the sum of two arbitrary columns. The expression will look like this: list(agg1 = sum(col1),agg2=sum(col2))

sum_expr = "list(agg1 = sum(%1, na.rm=TRUE), agg2 = sum(%2, na.rm=TRUE))"## We build the expression as a string

sum_expr = gsub("%1", "minute", sum_expr)
sum_expr = gsub("%2", "goal_number", sum_expr)
expr_sum = parse(text = sum_expr)## Parsing a string into an expression
GOALS[, eval(expr_sum)]
GOALS[, eval(expr_sum), by = penal]

#######
## We walked you through a few examples of datatable usage. This tutorial is not thorough but it should cover most of the basic usage.

## If you would like to know more about data.table implementation you can either set the "verbose" argument to TRUE

GOALS[, eval(expr_sum), by = penal, verbose = TRUE]

## Or set the "datatable.verbose" option to true:
## options(datatable.verbose=TRUE)
