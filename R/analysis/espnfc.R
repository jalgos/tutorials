
#Function that will read the csv file designated by fn and load it into an R datatable
#This is not the canonical way. The files are exports from mongo and all the columns have a character type. Conversion needs to be done
load_espnfc_table <- function(fn, cclasses)
{
    cat("loading file:", fn, "\n")
    D = as.data.table(read.table(fn, sep = ",", header = TRUE, colClasses = "character", quote = "\"", na.strings = c("nan.0", "-2147483648"), comment.char = "", encoding = "UTF-8"))
    names(cclasses) = names(D)
    convert_columns(D, cclasses)
    D
}

#Will plot a histogram and an empirical distribution of the goals by minute
goal_distribution <- function(F, gtype = "goals")
{
    if(missing(F)) F =GOALS[, penal != "penmiss" & minute < 90 & minute != 45 & is.na(extra_time)] #We need to filter out extra_time otherwise it gets confusing
    yl = "number of %1"
    title1 = "frequency of %1 per minute"
    title2 = "frequency of %1 per minute (1 bar per minute)"
    yl = gsub("%1", gtype, yl)
    title1 = gsub("%1", gtype, title1)
    title2 = gsub("%1", gtype, title2)
    dev.new()
    cat("plotting graph:", title1, "\n")
    GOALS[F, hist(minute,
                  breaks = 30,
                  ylab = yl,
                  main = title1)]#classic histogram
    press_key_to_continue()
    dev.new()
    cat("plotting graph:", title2, "\n")
    GOALS[F, hist(minute,
                  breaks = 90,
                  ylab = yl,
                  main = title2)]#classic histogram
    press_key_to_continue()
}

#Same same but with ggplot2
ggp2.goal_distribution <- function(F, gtype = "goals")
{
    if(missing(F)) F =GOALS[, penal != "penmiss" &  minute<90 & minute!=45 & is.na(extra_time)] #We need to filter out extra_time otherwise it gets confusing
    yl = "number of %1"
    title1 = "frequency of %1 per minute"
    title2 = "density of %1 per minute"
    yl = gsub("%1", gtype, yl)
    title1 = gsub("%1", gtype, title1)
    title2 = gsub("%1", gtype, title2)   
    dev.new()
    (ggplot(GOALS[F], aes(x = minute)) + geom_histogram() + ggtitle(title1) + ylab(yl))#I don't understand why this one doesn't plot
    dev.new()
    (ggplot(GOALS[F], aes(x = minute)) + geom_density() + ggtitle(title2))
}

#Function to filter valid boxscore. The criteria will be at least 1 shot in the game and the number of shots is greater than the number of goals
valid_bxs <- function()
{
    vgids = BXS[, list(goals = sum(g), shots = sum(sh)), by = list(game_id, team_id)][goals <= shots][shots != 0, unique(game_id)]
    BXS[, game_id %in% vgids]
}

bxs_analysis <- function()
{
    F = valid_bxs()
    #shot on target %, goal % 
    STP = BXS[F, list(shots = sum(sh), on_goal = sum(sg), goals = sum(g)), keyby = player_id]
    STP[PL, full_name := full_name]
    #who shoots on target:
    STP[, c("tgt_pct", "goal_pct") := list(on_goal / shots, goals / shots)]
    cat("Best 25 shooter on target\n")
    print(STP[shots > 100][order(tgt_pct, decreasing = TRUE)][1:25])
    press_key_to_continue()
    cat("Best 25 scorer (goal/shots ratio)\n")
    print(STP[shots > 100][order(goal_pct, decreasing = TRUE)][1:25])
    press_key_to_continue()
    #Goalies
    #A bit tricky here. We need to get the total number of shots on target achieved by the opponents
    setkey(BXS, game_id, team_id)
    setkey(GS, game_id, team1)
    BXS[GS, c("opp_sht", "saves") := list(shots_on_target2, shots_on_target2 - score2)]
    setkey(GS, game_id, team2)
    BXS[GS, c("opp_sht", "saves") := list(shots_on_target1, shots_on_target1 - score1)]
    SVS = BXS[F & opp_sht >= sv & pos == "G", list(saves = sum(saves), opp_shots = sum(opp_sht)), keyby = player_id]
    SVS[PL, c("full_name", "position") := list(full_name, as.character(position))]
    SVS[, save_pct := saves / opp_shots]
    cat("Best 25 goalkeeper (save percentage)\n")
    print(SVS[order(save_pct, decreasing = TRUE)][saves > 100][1:25])
    press_key_to_continue()
    cat("Hugo Lloris save percentage\n")
    print(SVS[full_name == "Hugo Lloris"])
    press_key_to_continue()
}

