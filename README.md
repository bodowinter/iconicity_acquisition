Statistical analysis for "Iconicity in the speech of children and adults"
=============

**Authors:** Lynn Perry, Marcus Perlman, Bodo Winter, Dominic Massaro, Gary Lupyan

## Libraries required for this analysis:

-	dplyr
-	readr
-	stringr
-	reshape2
-	car
-	mgcv
-	lme4
-	mgcv
-	itsadug

## Script files contained in this analysis:

-	analysis.R
	Contains all analyses
-	predict.glmm.R
	Function for retrieving 95% confidence intervals from mixed models

## Data files contained in this analysis:

-	iconicity.csv
	The main dataframe that contains iconicity ratings by word, concreteness etc.
-	word_selection.txt
	The subset of words analyzed in this study. Excludes words from Perry et al. (2015) and compounds.
-	onomatopoetic_words.txt
	List of onomatopoetic words excluded for some analyses.
-	ANC-spoken-lemma.txt
	American National Corpus frequency data
-	monaghan2014_systematicity.csv
	Monaghan et al. (2014)'s systematicity data.
-	childesPlusIconNEW.csv
	By-child CHILDES data (longitudinal)

