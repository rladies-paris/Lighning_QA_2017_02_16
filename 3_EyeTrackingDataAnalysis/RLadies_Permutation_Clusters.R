### SCRIPT BEGIN ##########################################################
# 15.2.2017 Katie Von Holzen
# This script is a quick introduction to 
# permutation clusters analysis using the tutorial vignettes
# from the eyetrackingR package
# also known as: Estimating Divergences
# http://www.eyetracking-r.com/vignettes/divergence

# suggested to first work through the 
# Preparing Your Data vignette
# http://www.eyetracking-r.com/vignettes/preparing_your_data
#
# and to read:
# Maris, E., Oostenveld, R., (2007). 
# Nonparametric statistical testing of EEG- 
# and MEG-data. Journal of Neuroscience 
# Methods 164 (1), 177–190.
# doi: 10.1016/j.jneumeth.2007.03.024

# and:
# Von Holzen, K., Mani, M., (2012).
# Language nonselective lexical access in 
# bilingual toddlers. Journal of Experimental 
# Child Psychology, 113, 569–586. 
# doi: 10.1016/j.jecp.2011.02.002
#
#
# On each trial, infants were shown a picture of an animate object 
# (e.g., a horse) and an inanimate object (e.g., a spoon). 
# After inspecting the images, they disappeared and they heard a 
# label referring to one of them (e.g., “The horse is nearby!”). 
# Finally, the objects re-appeared on the screen and they were 
# prompted to look at the target (e.g., “Look at the horse!”).
#
# look here at a time course plot

rm(list=ls())
library(eyetrackingR)
library(reshape)
library(ggplot2)
library(grid)
library(gridExtra)


# SET YOUR WORKING DIRECTORY NOW
if(Sys.info()['sysname'] == "Darwin") {
  setwd("~/Documents/R_Git/eyetrackingWorkshop")
  core <- "~/Documents/R_Git/eyetrackingWorkshop" 
  git_loc <-"~/Documents/R_Git/eyetrackingWorkshop"
} else {
  setwd("C:/Users/katie/R Git/eyetrackingWorkshop")
  core <- "C:/Users/katie/R Git/eyetrackingWorkshop"
  git_loc <- "C:/Users/katie/R Git/eyetrackingWorkshop"
}


############################################
# load data
############################################

data("word_recognition")
data <- make_eyetrackingr_data(word_recognition, 
                               participant_column = "ParticipantName",
                               trial_column = "Trial",
                               time_column = "TimeFromTrialOnset",
                               trackloss_column = "TrackLoss",
                               aoi_columns = c('Animate','Inanimate'),
                               treat_non_aoi_looks_as_missing = TRUE
)

# subset to response window post word-onset
response_window <- subset_by_window(data, 
                                    window_start_time = 15500, 
                                    window_end_time = 21000, 
                                    rezero = FALSE)

############################################
# clean data
############################################

# analyze amount of trackloss by subjects and trials
(trackloss <- trackloss_analysis(data = response_window))

# remove trials with > 25% of trackloss
response_window_clean <- clean_by_trackloss(data = response_window,
                                            trial_prop_thresh = .25)

response_window_clean$Target <- as.factor( ifelse(test = grepl('(Spoon|Bottle)', 
                                                  response_window_clean$Trial), 
                                                  yes = 'Inanimate', 
                                                  no  = 'Animate') )

############################################
# create and examine time course
############################################
response_time <- make_time_sequence_data(response_window_clean,
                                         time_bin_size = 100, 
                                         predictor_columns = c("Target"),
                                         aois = "Animate",
                                         summarize_by = "ParticipantName" )

# visualize timecourse
plot(response_time, predictor_column = "Target") + 
  theme_light() +
  coord_cartesian(ylim = c(0,1))+
  geom_hline(yintercept=0.5)

############################################
# analyze time course
############################################
tb_analysis <- analyze_time_bins(data = response_time, 
                                 predictor_column = "Target", 
                                 test = "t.test", alpha = .05)

plot(tb_analysis, type = "statistic") + theme_light()
# 55 T-TESTS!

############
# but, considering how many t-tests we run,
# we raise our chances of one of those tests 
# reaching our critical alpha

# what are our odds of getting a false alarm?

alpha <- .05
#alpha <- 0.0001666667
num_time_bins <- nrow(tb_analysis)
(prob_no_false_alarm_per_bin <- 1-alpha)

(prob_no_false_alarm_any_bin <- prob_no_false_alarm_per_bin^num_time_bins)

# our probability of at least one false alarm
(prob_at_least_one_false_alarm <- 1-prob_no_false_alarm_any_bin)

# 94% probability that at least one of these t-tests is a false alarm
# we are not controlling for the the family-wise error rate

################################################################################
################################################################################
# Bootstrapped cluster-based permutation analysis
################################################################################
################################################################################

# This analysis takes a summed statistic for each cluster of time bins 
# that pass some level of significance, and compares each to the "null" 
# distribution of sum statistics (obtained by bootstrap resampling data 
# within the largest of the clusters).
#
#######################################################################################################
#
# Advantages:
#
# 1. It naturally controls the false-alarm rate while sacrificing little sensitivity.
#
# 2. The implementation in eyetrackingR allows you to use this method with a variety 
# of statistical techniques (t.test, wilcox.test, lm, and lmer), so that continuous 
# predictors, covariates, etc. can also be included in the model being tested.
#
# Alternative:
#
# Bonferroni correction
# alpha/C (where C is the number of time samples)
# ex) 0.05/300
# Bonferroni corrected alpha: 0.0001666667
#
# but, neuroimaging tests are not independent from one another
# we do not necessarily know the spatial correlation of the data
# 
#######################################################################################################
# what does it do?
#
# 1. Runs a statistical test on each time-bin of your data.
#
# 2. Take the time-bins whose test passed the a threshold 
#    statistic (e.g., t > 2.26), and group them by adjacency. 
#    We will call these time-clusters.
#
# 3. For each time-cluster, calculate the sum of the statistics 
#    for the time-bins inside it.
#
# 4. Take the data inside the largest of these clusters. Randomly shuffle it.
#
# 5. Recalculate the sum-statistic for each of the shuffled datasets you sampled in (4).
#
# 6. Repeat steps (4) and (5) hundreds of times. This will lead to a distribution 
#    of summed-statistics, each representing the results of a statistical test on 
#    shuffled data. Intuitively, this distribution represents what kind of sum-statistics 
#    we would expect to get in a cluster by chance, if no effect were present 
#    (i.e., if the data were just randomly shuffled).
#
# 7. Compare the cluster-statistics from step (3) to the distribution found in (6) 
#    to obtain a p-value. So, for example, imagine we get a distribution of 
#    sum-statistics, and 6.4% of the sum-statistics in this distribution are 
#    greater than the sum-statistic of our cluster, then our p-value for this cluster is p = .064.
#######################################################################################################
#
# this test assumes our ignorance of the spatiotemporal locus of the effect
#
# The permutation test controls the FA rate in the following
# conditional sense: given the unordered set {D}={d}, under 
# exchangeability, the probability of observing a p-value that is
# less than the critical alpha-level is exactly equal to the critical alpha-level.
#
# so, we do comparisons multiple times, to determine how often we do get false
# positives with our data, and we then check to make sure that with our actual
# data we are above some critical alpha threshold (usually 0.05)`
#
# we need to set a threshold statistic

alpha = .05
num_sub = length(unique((response_window_clean$ParticipantName)))
threshold_t = qt(p = 1 - alpha/2, 
                 df = num_sub-1) # pick threshold t based on alpha = .05 two tailed

# now, we look for the clusters in the original data
# this is steps 1 and 2

df_timeclust <- make_time_cluster_data(response_time, 
                                       test= "t.test",
                                       paired = TRUE,
                                       predictor_column = "Target", 
                                       threshold = threshold_t) 

# so, where are the divergences in the data?
plot(df_timeclust) +
  ylab("T-Statistic")

summary(df_timeclust)
# between 16100 - 19300, cluster 1
# between 19400 - 20800, cluster 2

# step 3
# sum statistic
# so, the sum of the t-values within each cluster

# next is steps 4-6
# we are looking at the time clusters, provided by "df_timeclust"
# we are shuffling it up and then for each shuffle, recalculating the sum statistic for the time clusters
# we repeat this many times, which gives us a distribution of summed statistics
# if we had no effect in the data, we would expect 
#
# we want to see number of the sum-statistics in this distribution that are 
# greater than the sum-statistic of our cluster
#
# in this step, within our previously defined clusters, we shuffle the assignment of trials
# to new conditions, within subjects

system.time(clust_analysis <- analyze_time_clusters(df_timeclust, 
                                                    within_subj = TRUE, 
                                                    paired=TRUE, 
                                                    samples=150)) # in practice, you should use a lot more
# 150 samples
#    user  system elapsed 
#    45.79    0.02   46.14 
# 30 ms per sample

clust_analysis

plot(clust_analysis)

# x axis is actual summed t-statistic
# conclude that our conditions are significantly different from one another

############################################
# 1000 shuffles
############################################

load("clust_analysis1000.RData")

# not run:
#system.time(clust_analysis1000 <- analyze_time_clusters(df_timeclust, within_subj = TRUE, paired=TRUE, 
#                                                    samples=1000)) # in practice, you should use a lot more
# 1000 samples
# user  system elapsed 
# 301.37    0.01  301.97
# 30 ms per sample

clust_analysis1000

plot(clust_analysis1000)

.011*1000# .011*1000

# Comparisons using the time course analysis revealed that 
# the animate and inanimate conditions significantly deviated 
# from each other between 16100 and 19300 ms, cluster t statistic = 132.30, 
# Monte Carlo p < .0001, and between 19400 and 20800 ms, 
# cluster t statistic= 42.31, Monte Carlo p < .05.
#
# none of the 1000 random partitions resulted in a 
# cluster-level statistic that is larger in absolute value
#
# if we had done the test an infinite number of times,
# it would be a permutation p-value and not a Monte Carlo p
# so, proportion of random partitions in which the observed 
# test statistic is larger than the value drawn from the 
# permutation distribution

############################################
# variation for each test run
############################################


system.time(clust_analysis.a <- analyze_time_clusters(df_timeclust, within_subj = TRUE, paired=TRUE, 
                                                      samples=150)) # in practice, you should use a lot more

system.time(clust_analysis.b <- analyze_time_clusters(df_timeclust, within_subj = TRUE, paired=TRUE, 
                                                      samples=150)) # in practice, you should use a lot more

clust_analysis.a
clust_analysis.b

plot(clust_analysis.a)
plot(clust_analysis.b)

############################################
# simulate data!!!
############################################

# set number of participants
n <- 16 # number of participants
c.trials <- 6 # items per condition
t.length <- 5000 # trial length
prf <- 0.6 # 	Their preference between the two AOIs in the 
            #   "high" condition, where 1 is 100 preference). 
            #   In the "low" condition, their preference between 
            #   the two AOIs is equal, so default is no effect of condition.
ns <- NULL # specifying start and end of time-window in which there was 
           # substantial trackloss during the trial

sim.data <- simulate_eyetrackingr_data(num_participants = n,
                           num_items_per_condition = c.trials, 
                           trial_length = t.length, 
                           pref = prf,
                           pref_window = c(1, t.length), 
                           noisy_window = ns)
# AOI1 is correct looks in High condition
# AOI2 is correct looks in Low condition

# condition tells us which AOI was correct
# so, if condition is LOW and AOI1 is TRUE
# then this a look to the distractor
# if condition is LOW and AOI1 is FALSE
# this is a look to the target

# recode to have looks to Target be TRUE, depending on Condition
#sim.data$Target <- with(sim.data, ifelse(Condition == "Low" & AOI1 == "TRUE", FALSE,
#                               ifelse(Condition == "Low" & AOI1 == "FALSE", TRUE,
#                               ifelse(Condition == "High" & AOI1 == "TRUE", TRUE, FALSE))))


############################################
# create and examine time course
############################################
response_time.sim <- make_time_sequence_data(sim.data,
                                         time_bin_size = 100, 
                                         predictor_columns = c("Condition"),
                                         aois = "AOI1",
                                         #aois = "Target",
                                         summarize_by = "Participant" )

# visualize timecourse
plot(response_time.sim, predictor_column = "Condition") + 
  theme_light() +
  coord_cartesian(ylim = c(0,1))

############################################
# examine clusters
############################################

# we need to set a threshold statistic
alpha = .05
num_sub = n
threshold_t = qt(p = 1 - alpha/2, 
                 df = num_sub-1) # pick threshold t based on alpha = .05 two tailed

# now, we look for the clusters in the original data
# this is steps 1 and 2

df_timeclust <- make_time_cluster_data(response_time.sim, 
                                       test= "t.test",
                                       #paired = TRUE,
                                       predictor_column = "Condition", 
                                       threshold = threshold_t) 

# so, where are the divergences in the data?
plot(df_timeclust) +
  ylab("T-Statistic")

summary(df_timeclust)


# in this step, within our previously defined clusters, we shuffle the assignment of trials
# to new conditions, within subjects
system.time(clust_analysis <- analyze_time_clusters(df_timeclust, 
                                                    within_subj = FALSE, 
                                                    #paired=TRUE, 
                                                    samples=150)) # in practice, you should use a lot more
# examine cluster analysis
clust_analysis

# plot t values
plot(clust_analysis)

######
pref.08 <- clust_analysis
pref.06 <- clust_analysis
