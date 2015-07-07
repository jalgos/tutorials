#This example script will demonstrate how to use R to handle data
#It uses two very popular packages data.table and ggplot2
#It will output a few graphs and a few results using the ESPNFC dataset
#You should set up your directory structure as follow:
#two folders at the root of your project:
#-root/data: containing all the files present in the drive:https://drive.google.com/a/jalgos.com/?tab=mo#folders/0B1ej3OAh42ovbjVVTENQTnBWMzA

#-root/analysis: a folder analysis containing the present script and espnfc.R



source("analysis/util.R")
load_or_install("data.table")
load_or_install("ggplot2")
source("analysis/espnfc.R")
#The players table contains generic information about the players. This table can be joined against other table thanks to the player_id column
PL = load_espnfc_table("data/players.csv",c("character","numeric","factor","numeric","character","factor","character","factor","character"))
#indexing by player_id
setkey(PL,player_id)


#A mapping between team names and team_id.
TEAMS = load_espnfc_table("data/teams.csv",c("character","character","character"))
#indexing by team_id
setkey(TEAMS,team_id)


#Info about leagues. Mapping is done through league_id
LEAGUES = load_espnfc_table("data/leagues.csv",c("character","character"))
#indexing by league_id
setkey(LEAGUES,league_id)

#An entry per date with the number of games played and the number of leagues
DAYS = load_espnfc_table("data/days.csv",c("numeric","numeric","character"))

#A summary for each game played with basic statistics. Each record has a unique game_id
GS = load_espnfc_table("data/games.csv",c("factor","factor","factor","factor","numeric","factor","factor",rep("numeric",13),"factor",rep("numeric",8)))
#indexing by game_id
setkey(GS,game_id)

#A table containing all the goals, with the scorer the minute, the type (penalty, own goal etc...) if it was scored in extra time. Can be linked to other table through team_id, game_id, player_id
GOALS = load_espnfc_table("data/goals.csv",c("factor","factor","factor","numeric","integer","factor","factor","numeric"))

#fixing column penal
GOALS[,penal:=factor(gsub("\\s*$","",penal))]

#Table containing "advanced" aggregated stats. Shots on target, saves etc... One record per game per player.
BXS = load_espnfc_table("data/boxscore.csv",c("numeric","numeric","factor","factor","factor","numeric","factor",rep("numeric",7),"factor",rep("numeric",7)))

#All these variables are data tables and they are all global. It means they are accessible in any function

cat("Check out the data table tutorial in analysis/data_table_tutorial.R\n")

####Global consideration on goals

goal_distribution()

goal_distribution(GOALS[,penal =="pen" & minute!=45 &minute<90],"penalties")#distribution of penalties

goal_distribution(GOALS[,penal =="og" & minute!=45 &minute<90],"own goals")#distribution of own goals

goal_distribution(GOALS[,scorer =="Lionel Messi" & minute!=45 &minute<90],"Messi's goals")#distribution of Messi's goals


#25 best scorers of all time
NG = GOALS[,list(num_goals =.N),by = scorer]
cat("Top 25 scorers of all time\n")
print(NG[order(num_goals,decreasing=TRUE)][1:25])
press_key_to_continue()
   

bxs_analysis()
