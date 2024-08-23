
{smcl}

{* *! version 1.0.0 16Jun2021}
{cmd:help pmcstat}
{hline}

{title:Title}

{phang}
{bf:pmcstat} {hline 2} Calculates the C-statistic for a prediction model

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:pmcstat} {varlist},  
	[{opt han:ley}
	{opt mat:rix(name)}
	{opt noprint}]

{pstd}
where {varlist} must consist of two variables; the first variable must contain the linear predictor values from the model. 
The second variable must represent the event indicator i.e. the observed outcome (a binary variable coded 0 for non-events and 1 for events).
{p_end}

{synoptset 30}{...}
{synopthdr:pmcstat_options}
{synoptline}
{synopt :{opt han:ley}}use SE formula of Hanley (default uses Newcombe formula){p_end}
{synopt :{opt mat:rix(name)}}specifies a name for the output results matrix{p_end}
{synopt :{opt noprint}}suppress output{p_end}

{synoptline}

{marker description}{...}
{title:Description}

{pstd}{cmd: pmcstat} computes the C-statistic for a clinical prediction model at either model development or evaluation.
{cmd: pmcstat} can currently be used to calculate the C-statistic for models with binary outcomes. Survival (time-to-event) outcomes will be available in future versions. 
The aim of the package is faster computation of the C-statistic, which can take quite some time using alternative packages in Stata.
For example, {cmd: pmcstat} is ~140 times faster than {cmd: roctab} when computing the C-statistic on a sample of ~20k.
Latest version can be accessed from GitHub using: {cmd: }{stata "net from https://joieensor.github.io/pm-suite/"}

{marker options}{...}
{title:Options}

{phang}{opt hanley} specifies that the Hanley formula for standard error calculation be used.
Default calculation of the standard error uses the Newcombe formula. 

{phang}{opt matrix(name)} specifies a name for the results matrix to be stored.

{phang}{opt noprint} suppress the display of results output table.


{marker examples}{...}
{title:Examples}

{phang}
Below we will fit a simple prediction model to predict low birthweight.
For illustration only we will split our sample into train/test sets, while this is not recommended, 
it will allow us to illustrate the use of {cmd: pmcstat} for both development and validation scenarios. 

{pstd}Load dataset{p_end}
{phang2}{cmd:. }{stata webuse lbw, clear}{p_end}

{pstd}For illustration only, create a binary variable to define a development and validation cohort{p_end}
{phang2}{cmd:. }{stata gen val = runiform()<.4}{p_end}

{pstd}Derive a prediction model for low birth weight in the development cohort{p_end}
{phang2}{cmd:. }{stata logistic low age lwt i.race smoke ptl ht ui if val==0}{p_end}

{pstd}Generate a new variable containing the linear predictor values from the model for all individuals{p_end}
{phang2}{cmd:. }{stata predict lp, xb}{p_end}

{pstd}Use {cmd: pmcstat} to calculate the apparent C-statistic performance of the development model in the development cohort{p_end}
{phang2}{cmd:. }{stata pmcstat lp low if val==0}{p_end}

{pstd}Now calculate the models C-statistic performance in the validation cohort{p_end}
{phang2}{cmd:. }{stata pmcstat lp low if val==1}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:pmcstat} stores the following in {cmd:r()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(cstat)}}C-statistic for model discrimination{p_end}
{synopt:{cmd:r(cstat_se)}}Standard error for C-statistic{p_end}
{synopt:{cmd:r(cstat_lb)}}Lower bound of 95% CI for C-statistic{p_end}
{synopt:{cmd:r(cstat_ub)}}Upper bound of 95% CI for C-statistic{p_end}
{synopt:{cmd:r(cord)}}Number of concordant pairs in calculations{p_end}
{synopt:{cmd:r(disc)}}Number of discordant pairs in calculations{p_end}
{synopt:{cmd:r(ties)}}Number of ties in calculations{p_end}
{synopt:{cmd:r(obs)}}Total observations used in calculations{p_end}
{p2colreset}{...}

{title:Authors}

{phang}Joie Ensor, University of Birmingham {break}
j.ensor@bham.ac.uk{p_end}

{title:Acknowledgements}

{phang}With thanks to Emma Martin{p_end}


{title:Also see}

{psee}
Online: {helpb pmsampsize}, {helpb pmvalsampsize}, {helpb pmcalplot}
{p_end}
