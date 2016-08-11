## Reboot of the iconicity ~ AOA analysis
## August 5, 2016
## Edit August 8-12, 2016
## Authors: Lynn K. Perry, Marcus Perlman, Bodo Winter, ...
## ... Dominic W. Massaro, Gary Lupyan


##------------------------------------------------------------------
## Pre-processing:
##------------------------------------------------------------------

## Libraries:

library(car)
library(lme4)
library(mgcv)
library(readr)
library(stringr)
library(reshape2)
library(dplyr)

## Set global options:

options(stringsAsFactors = F)
options(dplyr.width = Inf)

## Load in data:

setwd('/Users/teeniematlock/Desktop/research/iconicity/AOA_paper/analysis/data/')
icon <- read.csv('iconicity.csv') %>% tbl_df()
# 'down' was redundantly collected in both list1 and list2
# here, the average rating for both 'down' was used

## Get only those words that Lynn wants to analyze:
# (this excludes words from list1, as well as several compounds etc.)

these_words <- dget('word_selection.txt')
icon <- icon %>% filter(Word %in% these_words)

## Rename and log-transform frequencies:

icon <- icon %>% rename(POS = SUBTLEX_dom_POS,
	Freq = SUBTLEX_Rawfreq) %>%
		mutate(LogFreq = log10(Freq))

## Retrieve words with missing POS tags:

POS <- icon %>%
	filter(is.na(POS) | POS %in% c('#N/A', 'Unclassified')) %>%
		select(Word)
POS$POS <- c('Verb', 'Verb', 'Verb', 'Verb',
	'Interjection', 'Noun', 'Onomatopoetic',
	'Onomatopoetic', 'Onomatopoetic', 'Onomatopoetic',
	'Onomatopoetic', 'Noun', 'Adjective', 'Adjective')

## Merge missing words into there:

icon[match(POS$Word, icon$Word), ]$POS <- POS$POS

## Lump together lexical categories into meaningful subcategories:

gram <- c('Article', 'Conjunction', 'Determiner',
	'Ex', 'Not', 'Number', 'Pronoun', 'Preposition', 'To')
icon[icon$POS %in% gram, ]$POS <- 'Grammatical'

## Fix the 'Name' category:
# (includes a lot of things that are definitely not names)

nam <- icon %>% filter(POS == 'Name') %>% select(Word)
nam$POS <- c('Adjective', 'Noun', 'Noun', 'Interjection',
	'Noun', 'Verb', 'Adjective', 'Noun', 'Noun',
	'Noun', 'Adjective', 'Noun', 'Verb', 'Adjective',
	'Adjective', rep('Noun', 17), 'Adjective', 'Adjective',
	rep('Noun', 8), 'Verb', 'Noun')
icon[icon$POS %in% 'Name', ]$POS <- nam$POS

## Merge onomatopoetic and interjections together:

icon[icon$POS == 'Onomatopoetic', ]$POS <- 'Interjection'

## Add onomatopoeia information to dataset:

onom <- dget('onomatopoetic_words.txt')
icon$Onomatopoetic <- 0
icon[icon$Word %in% onom, ]$Onomatopoetic <- 1

## Merge American National Corpus frequency in there:

ANC <- read.table('ANC-spoken-lemma.txt',
	header = F, skip = 13) %>% tbl_df()
ANC <- ANC %>% rename(Word = V1, Lemma = V2,
	POS = V3, Freq = V4)

## Take aggregates by word form:

ANC <- ANC %>% group_by(Word) %>% summarise(Freq = sum(Freq))

## Log-transform this:

ANC <- ANC %>% mutate(LogFreq = log10(Freq))

## Merge into iconicity data frame:

icon$ANC_LogFreq <- ANC[match(icon$Word, ANC$Word), ]$LogFreq



##------------------------------------------------------------------
## Analysis of age of acquisition ratings:
##------------------------------------------------------------------

## Does iconicity predict age of acquisition?

xmdl.AOA <- lm(KupermanAOA ~ Iconicity + Conc +
	POS + NMorph + ANC_LogFreq,
	data = icon)
summary(xmdl.AOA)	# yes!
anova(xmdl.AOA)		# yes!

## Model validation:

plot(fitted(xmdl.AOA),
	residuals(xmdl.AOA))	# o.k.; slight heteroskedasticity
qqnorm(residuals(xmdl.AOA))
qqline(residuals(xmdl.AOA))	# slight positive skew

## Check for collinearity:

vif(xmdl.AOA)	# no problem

## Exclude onomatopoetic words:

noOnom <- icon[!as.logical(icon$Onomatopoetic), ]
nrow(icon) - nrow(noOnom)	# excluded 27 datapoints

## Is there iconicity ~ AOA even without onomatopoetics?

xmdl.noOnom <- lm(KupermanAOA ~ Iconicity + Conc +
	POS + NMorph + ANC_LogFreq,
	data = noOnom)
summary(xmdl.noOnom)	# yes!
anova(xmdl.noOnom)		# yes!

## Incorporate systematicity from Monaghan et al. (2014):

sys <- read.csv('monaghan2014_systematicity.csv')
icon$Sys <- sys[match(icon$Word, sys$word), ]$relativeIconicity

## Is there iconicity ~ AOA controlling for systematicity?

xmdl.sys <- lm(KupermanAOA ~ Iconicity + Conc +
	POS + NMorph + ANC_LogFreq + Sys,
	data = icon)
summary(xmdl.sys)	# yes!
anova(xmdl.sys)		# yes!

## This also all works with SUBTLEX frquencies!



##------------------------------------------------------------------
## Added variable plot of iconicity after residualizing other effects:
##------------------------------------------------------------------

## Model without iconicity:

all_compl <- with(icon, complete.cases(KupermanAOA, Conc, POS,
	NMorph, ANC_LogFreq, Iconicity))
icon.noNA <- icon[all_compl, ]
xmdl.noicon <- lm(KupermanAOA ~ Conc +
	POS + NMorph + ANC_LogFreq,
	data = icon.noNA)
xmdl <- lm(residuals(xmdl.noicon) ~ Iconicity, data = icon.noNA)

## Extract predictions:

xvals <- seq(from = -7, to = 7, by = 0.001)
xpred <- data.frame(Iconicity = xvals)
xpred <- as.data.frame(predict(xmdl, newdata = xpred, se.fit = T)[1:2])
xpred$UB <- xpred$fit + 1.96 * xpred$se.fit
xpred$LB <- xpred$fit - 1.96 * xpred$se.fit

## Vector of words to display:

these_words <- c('beep', 'click', 'moo', 'icky', 'roar', 'heaven', 'scratchy',
	'pajamas', 'silent', 'would', 'hamster', 'peekaboo',
	'quality', 'hum', 'drag', 'staff', 'incentive', 'computer',
	'scale', 'drag', 'buzz', 'pretend', 'clamp', 'canoe', 'mushy',
	'socialist', 'bureau')
these_words <- icon.noNA$Word %in% these_words

## Make a plot of this:

quartz('', 9, 6.5)
par(mai = c(1.5, 1.5, 0.25, 0.25))
plot(1, 1, type = 'n',
	xlim = c(-4, 6), ylim = c(-6, 8),
	xlab = '', ylab = '',
	xaxt = 'n', yaxt = 'n')
points(x = icon.noNA[!these_words, ]$Iconicity,
	y = residuals(xmdl.noicon)[!these_words],
	pch = 21, cex = 1.25,
	bg = rgb(0, 0, 0, 0.4), col = NA)
text(x = icon.noNA[these_words, ]$Iconicity,
	y = residuals(xmdl.noicon)[these_words],
	labels = icon.noNA[these_words, ]$Word,
	pch = 21, cex = 1.1,
	bg = rgb(0, 0, 0, 0.4), col = NA)
polygon(x = c(xvals, rev(xvals)),
	y = c(xpred$UB, rev(xpred$LB)),
	border = NA, col = rgb(0, 0, 0, 0.4))
abline(xmdl, lwd = 2)
## Add axes:
axis(side = 1, at = seq(-4, 6, 2), font = 2,
	lwd.ticks = 2, cex.axis = 1.25)
mtext(text = 'Average Iconicity Rating', side = 1,
	line = 4, font = 2, cex = 2)
axis(side = 2, at = seq(-6, 8, 2), font = 2,
	lwd.ticks = 2, cex.axis = 1.25, las = 2)
mtext(text = 'Age of Acquisition Rating', side = 2,
	line = 5, font = 2, cex = 2)
mtext(text = '(Residuals)', side = 2, line = 3.15,
	font = 2, cex = 1.5)
box(lwd = 2)




##------------------------------------------------------------------
## Analysis of child & adult production frequencies:
##------------------------------------------------------------------

## Get subset of relevant variables:

childes <- icon %>% filter(!is.na(ChildesAge12Freq)) %>% 
	select(Word, Iconicity, ParentalInputRawFreq,
	Onomatopoetic,
	LogFreq, ANC_LogFreq,
	Conc, NMorph, POS,
	ChildesAge12Freq:ChildesAge69Freq)

## Reshape into long format:

childes <- melt(childes, id.vars = c('Word', 'Iconicity', 'ParentalInputRawFreq',
	'Onomatopoetic', 'LogFreq', 'ANC_LogFreq', 'Conc', 'NMorph', 'POS')) %>% tbl_df()

## Rename:

childes <- childes %>% rename(Age = variable, Freq = value)

## Log-transform parental input frequency and children's frequency:

childes <- childes %>% mutate(LogParent = log10(ParentalInputRawFreq),
	LogChild = log10(Freq + 1))

## Clean up the age column:

childes <- childes %>%
	mutate(Age = as.numeric(str_extract(as.character(Age), '[0-9]+')))

## Since we are going to fit an interaction, we should center:

childes <- childes %>% mutate(Age_c = Age - mean(Age),
	Iconicity_c = Iconicity - mean(Iconicity),
	Age_z = Age_c / sd(Age_c),
	Iconicity_z = Iconicity / sd(Iconicity),
	Conc_c = Conc - mean(Conc, na.rm = T),
	NMorph_c = NMorph - mean(NMorph, na.rm = T),
	ANC_LogFreq_c = ANC_LogFreq - mean(ANC_LogFreq, na.rm = T),
	LogParent_c = LogParent - mean(LogParent, na.rm = T))

## Make a Poisson model out of this:

xmdl <- lmer(LogChild ~ Iconicity_z + Age_z + Iconicity_z:Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + POS + 
	(1 + Age_z|Word), data = childes, REML = F)
summary(xmdl)	# inspect

## Model validation:

plot(fitted(xmdl), residuals(xmdl))
qqnorm(residuals(xmdl))
qqline(residuals(xmdl))

## Check variance inflation facors in corresponding linear model:

childmdl <- lm(LogChild ~ Iconicity_z + Age_z + Iconicity_z:Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + POS,
	data = childes)
vif(childmdl)	# no reason for concern, all < 3

## Construct comparison models for likelihood ratio tests:

xmdl.noint <- lmer(LogChild ~ Iconicity_z + Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + POS + 
	(1 + Age_z|Word), data = childes, REML = F)
xmdl.noage <- lmer(LogChild ~ Iconicity_z + Iconicity_z:Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + POS + 
	(1 + Age_z|Word), data = childes, REML = F)
xmdl.noicon <- lmer(LogChild ~ Age_z + Iconicity_z:Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + POS + 
	(1 + Age_z|Word), data = childes, REML = F)

## Perform likelihood ratio tests:

anova(xmdl.noint, xmdl, test = 'Chisq')
anova(xmdl.noage, xmdl, test = 'Chisq')	# p < 0.001
anova(xmdl.noicon, xmdl, test = 'Chisq')	# p < 0.001

## Check this for whether onomatopoeitic words matter:

childes_noOnom <- childes %>% filter(Onomatopoetic == 0)
xmdl.noOnom <- lmer(LogChild ~ Iconicity_z + Age_z + Iconicity_z:Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + POS + 
	(1 + Age_z|Word), data = childes_noOnom, REML = F)
xmdl.noOnom.noint <- lmer(LogChild ~ Iconicity_z + Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + POS + 
	(1 + Age_z|Word), data = childes_noOnom, REML = F)
anova(xmdl.noOnom.noint, xmdl.noOnom, test = 'Chisq')

## Make a 3D plot of this:
## (based on code by Roger Mundry)

library(rsm)
mdl.red <- lm(LogChild ~ Age_z * Iconicity_z, data = childes)
persp(mdl.red, Age_z ~ Iconicity_z, zlab = 'LogChild',
	col = 'lightgrey', theta = 45, phi = 25)

## Fit a model of iconicity on adult frequencies:

icon <- icon %>% mutate(Conc_c = Conc - mean(Conc, na.rm = T),
	NMorph_c = NMorph - mean(NMorph, na.rm = T),
	Iconicity_z = (Iconicity - mean(Iconicity, na.rm = T)) / sd(Iconicity, na.rm = T))
summary(xmdl.adult <- lm(ANC_LogFreq ~ Iconicity_z + 
	Conc_c + NMorph_c + POS, data = icon))
anova(xmdl.adult)

## Fit a model of iconicity on parental input frequencies:

icon <- icon %>% mutate(LogParent = log10(ParentalInputRawFreq),
	ANC_LogFreq_c = ANC_LogFreq - mean(ANC_LogFreq, na.rm = T))
summary(xmdl.parent <- lm(LogParent ~ Iconicity_z + ANC_LogFreq_c +
	Conc_c + NMorph_c + POS, data = icon))
anova(xmdl.parent)
summary(xmdl.parent.nocontrol <- lm(LogParent ~ Iconicity_z + 
	Conc_c + NMorph_c + POS, data = icon))
anova(xmdl.parent.nocontrol)

## Model validation:

plot(fitted(xmdl.adult), residuals(xmdl.adult))
plot(fitted(xmdl.parent), residuals(xmdl.parent))
qqnorm(residuals(xmdl.adult))
qqline(residuals(xmdl.adult))	# perfect
qqnorm(residuals(xmdl.parent))
qqline(residuals(xmdl.parent))	# perfect




##------------------------------------------------------------------
## Create a plot of production frequencies:
##------------------------------------------------------------------

## For plotting adult & parental frequencies, models that ignore POS:

summary(plot.adult <- lm(ANC_LogFreq ~ Iconicity_z + 
	Conc_c + NMorph_c, data = icon))
summary(plot.parent <- lm(LogParent ~ Iconicity_z + ANC_LogFreq_c +
	Conc_c + NMorph_c, data = icon))

## Get predictions for that:

xvals <- seq(from = -2.5, to = 5.0, length.out = 1000)
xpred <- data.frame(Iconicity_z = xvals,
	Conc_c = 0, NMorph_c = 0, ANC_LogFreq_c = 0)
adult.pred <- as.data.frame(predict.lm(plot.adult,
	newdata = xpred, se.fit = T)[1:2])
parent.pred <- as.data.frame(predict.lm(plot.parent,
	newdata = xpred, se.fit = T)[1:2])
adult.pred <- adult.pred %>% mutate(UB = fit + 1.96 * se.fit,
	LB = fit - 1.96 * se.fit)
parent.pred <- parent.pred %>% mutate(UB = fit + 1.96 * se.fit,
	LB = fit - 1.96 * se.fit)

## For plotting CHILDS frequencies fit a model that ignores POS:

xmdl.plot <- lmer(LogChild ~ Iconicity_z + Age_z + Iconicity_z:Age_z + 
	Conc_c + NMorph_c + ANC_LogFreq_c + LogParent_c + 
	(1 + Age_z|Word), data = childes, REML = F)

## Find z-scores corresponding to age = 12, age = 39 and age = 69:

age_vals <- unique(childes[childes$Age %in% c(12, 39, 69), ]$Age_z)

## Create an empty data frame for predictions:

xpred <- data.frame(Iconicity_z = rep(xvals, 3),
	Age_z = rep(age_vals, each = length(xvals)),
	Conc_c = 0, NMorph_c = 0, ANC_LogFreq_c = 0, LogParent_c = 0)

## Source prediction function:

source('/Users/teeniematlock/Desktop/research/iconicity/AOA_paper/analysis/predict.glmm.R')

## Get predictions:

xpred <- predict.glmm(xmdl.plot, newdata = xpred)

## Restrict to observed range:

xpred <- xpred %>% filter(Iconicity_z >= min(childes$Iconicity_z),
	Iconicity_z <= max(childes$Iconicity_z))

age12 <- xpred[xpred$Age_z == age_vals[1], ]
age39 <- xpred[xpred$Age_z == age_vals[2], ]
age69 <- xpred[xpred$Age_z == age_vals[3], ]

## Likewise, restrict adult & parent to CHILDES range for plotting purposes:

mycondition <- xvals >= min(icon$Iconicity_z) & 
	xvals <= max(icon$Iconicity_z)
adult.pred <- adult.pred[mycondition, ]
parent.pred <- parent.pred[mycondition, ]
xvals_red <- xvals[mycondition]

## Make a plot:

quartz('', 11, 6)
par(mai = c(0, 0.25, 0.25, 0.25),
	omi = c(1.5, 1.5, 1, 1.5),
	mfrow = c(1, 2))
plot(1, 1, type = 'n',
	xlim = c(-2.5, 5), ylim = c(0, 1.25),
	xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
mtext('CHILDES Child Frequency', line = 1, font = 2, cex = 1.5)
text(labels = '(a)',
	x = -2.2, y = 1.25 - (0.03 * 1.25),
	font = 2, cex = 1.5)
# Plot lines, CHILDES:
points(x = age12$Iconicity_z[-c(1:2, nrow(age12) - 1, nrow(age12))],
	y = age12$LogChild[-c(1:2, nrow(age12) - 1, nrow(age12))],
	type = 'l', lwd = 4)
polygon(x = c(age12$Iconicity_z, rev(age12$Iconicity_z)),
	y = c(age12$UB, rev(age12$LB)),
	border = NA, col = rgb(0, 0, 0, 0.4))
points(x = age69$Iconicity_z[-c(1:2, nrow(age69) - 1, nrow(age69))],
	y = age69$LogChild[-c(1:2, nrow(age69) - 1, nrow(age69))],
	type = 'l', lwd = 4)
polygon(x = c(age69$Iconicity_z, rev(age69$Iconicity_z)),
	y = c(age69$UB, rev(age69$LB)),
	border = NA, col = rgb(0, 0, 0, 0.4))
## Add labels:
text(x = 2.25,
	y = 0.94,
	labels = 'Month = 12',
	font = 2, cex = 1.25, adj = c(0.5, NA),
	srt = 26)
text(x = 2.25,
	y = 0.58,
	labels = 'Month = 69',
	font = 2, cex = 1.25, adj = c(0.5, NA),
	srt = -9)
# Axes:
axis(side = 1, at = seq(-2.5, 5, 2.5), font = 2,
	lwd.ticks = 2, cex.axis = 1.25)
axis(side = 2, at = seq(0, 1.25, 0.25), font = 2,
	lwd.ticks = 2, cex.axis = 1.25, las = 2)
mtext(text = 'Log Frequency', side = 2,
	line = 4.15, font = 2, cex = 2)
box(lwd = 2)
# Plot 2:
plot(1, 1, type = 'n',
	xlim = c(-2.5, 5), ylim = c(0.5, 2.5),
	xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
mtext('Adult Frequency', line = 1, font = 2, cex = 1.5)
text(labels = '(b)',
	x = -2.2, y = 2.5 - (0.03 * 2),
	font = 2, cex = 1.5)	
# Plot lines, adults:
points(x = xvals_red, y = adult.pred$fit,
	type = 'l', lwd = 4)
polygon(x = c(xvals_red, rev(xvals_red)), 
	y = c(adult.pred$UB, rev(adult.pred$LB)),
	border = NA, col = rgb(0, 0, 0, 0.4))
points(x = xvals_red, y = parent.pred$fit,
	type = 'l', lwd = 4)
polygon(x = c(xvals_red, rev(xvals_red)), 
	y = c(parent.pred$UB, rev(parent.pred$LB)),
	border = NA, col = rgb(0, 0, 0, 0.4))
## Add labels:
text(x = 1.35,
	y = 2.03,
	labels = 'Child-directed',
	font = 2, cex = 1.25, adj = c(0.5, NA),
	srt = 27)
text(x = 1.35,
	y = 1,
	labels = 'Adult-directed',
	font = 2, cex = 1.25, adj = c(0.5, NA),
	srt = -42)
# Axes:
axis(side = 1, at = seq(-2.5, 5, 2.5), font = 2,
	lwd.ticks = 2, cex.axis = 1.25)
axis(side = 4, at = seq(0.5, 2.5, 0.5), font = 2,
	lwd.ticks = 2, cex.axis = 1.25, las = 2)
box(lwd = 2)
title(main = 'Average Iconicity Rating\n(z scores)',
	outer = TRUE, line = -24, cex.main = 2, font.main = 2)




##------------------------------------------------------------------
## Analysis of CHILDES data longitudinally:
##------------------------------------------------------------------

## Load in CHILDES frequency:

setwd('/Users/teeniematlock/Desktop/research/iconicity/AOA_paper/analysis/data/')
CHILD <- read_csv('childesPlusIconNEW.csv') %>% tbl_df()

## Check counts:

length(unique(CHILD$name))	# 323 children
table(CHILD$name)

## Add covariates:

these_columns <- c('ANC_LogFreq_c', 'POS', 'NMorph_c', 'Conc_c', 'Onomatopoetic')
CHILD <- cbind(CHILD,
	childes[match(CHILD$Word, childes$Word), these_columns])

## Get rid of all children with less than 30 data points:

child_long <- table(CHILD$name)
not_these <- names(child_long[child_long <= 30])
CHILD <- CHILD %>% filter(!name %in% not_these)

## Get only those kids for which there is longitudinal data, at least five time points:

child_long <- table(CHILD$name, CHILD$age)
child_long <- apply(child_long, 1, FUN = function(x) sum(x != 0))
not_these <- names(child_long[child_long <= 5])
CHILD <- CHILD %>% filter(!name %in% not_these)

## Sum frequencies over sessions:

CHILD_freqs <- CHILD %>%
	group_by(Word, rating, ANC_LogFreq_c, POS, NMorph_c, Onomatopoetic, 
			Conc_c, name, age) %>%
		summarise(freq = sum(rawFreq))

## Append cumulative frequencies per age:

CHILD_freqs <- CHILD_freqs %>%
	group_by(name, ANC_LogFreq_c, POS, NMorph_c, Conc_c, Onomatopoetic,
		age) %>%
	summarise(freq_all = sum(freq)) %>% right_join(CHILD_freqs)

## Compute weighted means:

CHILD_means <- CHILD_freqs %>%
	group_by(name, age) %>%
		summarise(iconicity = weighted.mean(rating, freq),
			NMorph_c = weighted.mean(NMorph_c, freq, na.rm = T),
			Conc_c = weighted.mean(Conc_c, freq, na.rm = T),
			ANC_LogFreq_c = weighted.mean(ANC_LogFreq_c, freq, na.rm = T),
			Onomatopoetic = weighted.mean(Onomatopoetic, freq, na.rm = T))

## Make a line plot of this:

plot(1, 1, xlim = c(12, 117),
	ylim = c(-1, 3.5), type = 'n')
all_children <- unique(CHILD_means$name)
for (i in 1:length(all_children)) {
	this_child <- all_children[i]
	this_df <- filter(CHILD_means, name == this_child)
	points(this_df$age, this_df$iconicity,
		type = 'l', lwd = 2, col = rgb(0, 0, 0, 0.4))
	}

## Compute 100% of completion (this normalizing by child for different observation lengths):

CHILD_means <- CHILD_means %>%
	group_by(name) %>%
		mutate(age01 = (age - min(age)) / diff(range(age)))

## Make a line plot of this:

plot(1, 1, xlim = c(0, 1),
	ylim = c(-1, 3.5), type = 'n')
all_children <- unique(CHILD_means$name)
for (i in 1:length(all_children)) {
	this_child <- all_children[i]
	this_df <- filter(CHILD_means, name == this_child)
	points(this_df$age01, this_df$iconicity,
		type = 'l', lwd = 2, col = rgb(0, 0, 0, 0.4))
	}

## Make a simple mixed model of this:

CHILD_means <- CHILD_means %>%
	mutate(age_c = age - mean(age),
		age_z = age_c / sd(age_c))
summary(xmdl <- lmer(iconicity ~ age_z + (1 + age_z|name),
	data = CHILD_means, REML = F))
xmdl.null <- lmer(iconicity ~ 1 + (1 + age_z|name),
	data = CHILD_means, REML = F)
anova(xmdl.null, xmdl)

## Same for 0 to 1 range:

summary(xmdl <- lmer(iconicity ~ age01 + (1 + age01|name),
	data = CHILD_means, REML = F))
xmdl.null <- lmer(iconicity ~ 1 + (1 + age01|name),
	data = CHILD_means, REML = F)
anova(xmdl.null, xmdl)

## Check random effects:

sum(coef(xmdl)$name[, 2] > 0)	# 10 children had positive slopes
nrow(coef(xmdl)$name)	# of 150!

## With covariates:

summary(xmdl.covar <- lmer(iconicity ~ age01 +
	NMorph_c + ANC_LogFreq_c + Conc_c + 
	(1 + age01|name),
	data = CHILD_means, REML = F))
xmdl.covar.null <- lmer(iconicity ~ 1 +
	NMorph_c + ANC_LogFreq_c + Conc_c +
	(1 + age01|name), 
	data = CHILD_means, REML = F)
anova(xmdl.covar.null, xmdl.covar, test = 'Chisq')

## A model with no onomatopoetic words / interjections:

CHILD_means_onom <- CHILD_freqs %>%
	filter(Onomatopoetic == 0) %>%
	group_by(name, age) %>%
		summarise(iconicity = weighted.mean(rating, freq),
			NMorph_c = weighted.mean(NMorph_c, freq, na.rm = T),
			Conc_c = weighted.mean(Conc_c, freq, na.rm = T),
			ANC_LogFreq_c = weighted.mean(ANC_LogFreq_c, freq, na.rm = T))
CHILD_means_onom <- CHILD_means_onom %>%
	group_by(name) %>%
		mutate(age01 = (age - min(age)) / diff(range(age)))
summary(xmdl <- lmer(iconicity ~ age01 + (1 + age01|name),
	data = CHILD_means_onom, REML = F))
xmdl.null <- lmer(iconicity ~ 1 + (1 + age01|name),
	data = CHILD_means_onom, REML = F)
anova(xmdl.null, xmdl)

## A smaller model for those hat have at least 10 time points:

child_long <- table(CHILD_means$name, CHILD_means$age)
child_long <- apply(child_long, 1, FUN = function(x) sum(x != 0))
not_these <- names(child_long[child_long <= 10])
CHILD_means_red <- CHILD_means %>% filter(!name %in% not_these)

summary(xmdl.red <- lmer(iconicity ~ age01 + (1 + age01|name),
	data = CHILD_means_red, REML = F))
xmdl.null <- lmer(iconicity ~ 1 + (1 + age01|name),
	data = CHILD_means_red, REML = F)
anova(xmdl.null, xmdl.red)

## Check random effects:

sum(coef(xmdl.red)$name[, 2] > 0)	# 0 children had positive slopes
nrow(coef(xmdl.red)$name)	# of 66!

## Plot the random effects estimates:

quartz('', 9, 6.5)
par(mai = c(1.5, 1.75, 0.25, 0.25))
plot(1, 1, type = 'n',
	xlim = c(0, 1), ylim = c(0, 1.5),
	xlab = '', ylab = '',
	xaxt = 'n', yaxt = 'n')
all_children <- unique(CHILD_means_red$name)
for (i in 1:length(all_children)) {
	this_df <- unlist(coef(xmdl.red)$name[i, ])
	abline(a = this_df[1],
		b = this_df[2],
		lwd = 2, col = rgb(0, 0, 0, 0.4))
	}
## Add axes:
axis(side = 1, at = seq(0, 1, 0.25),
	font = 2, lwd.ticks = 2, cex.axis = 1.25)
mtext(text = 'Normalized Time', side = 1,
	line = 4, font = 2, cex = 2)
axis(side = 2, at = seq(0, 1.5, 0.25), font = 2,
	lwd.ticks = 2, cex.axis = 1.25, las = 2)
mtext(text = 'Average Iconicity', side = 2,
	line = 5.5, font = 2, cex = 2)
mtext(text = '(Frequency-weighted)', side = 2,
	line = 3.45, font = 2, cex = 1.5)
box(lwd = 2)

## For GAM analysis we need more time points
## to reliably assess nonlinearity:

library(mgcv)
CHILD_means_red$name <- as.factor(CHILD_means_red$name)
summary(xgam <- bam(iconicity ~ s(age01, k = 5) +
	s(age01, name, bs = 'fs', k = 5, m = 1),
	data = CHILD_means_red, method = 'fREML'))
# small edf of age01 suggests near-linear effect

## Plot GAM:

library(itsadug)
plot_smooth(xgam, 'age01', rm.ranef = T, main = 'Effect of t(ime)')


