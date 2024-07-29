
{smcl}

{* *! version 1.0.0 23Jul2024}
{cmd:help pmstabilityss}
{hline}

{title:Title}

{phang}
{bf:pmstabilityss} {hline 2} Calculates the minimum sample size required for developing a prediction model targeting precise individual risk estimates

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:pmstabilityss} {indepvars}, {opt prev:alence(real)} 
	[{opt cstat:istic(real 0)}
	{opt beta:s(numlist)}
	{opt n(numlist)}
	{opt st:andardised}  
	{opt lp(varname)} 
	{opt sub:group(varname)} 
	{opt simobs(int 500000)}
	{opt pci:width(numlist)}
	{opt pcut:points(numlist ascending >=0 <=1)}
	{opt thresh:old(real 0)}
	{opt tol:erance(real 0.005)}
	{opt logitpvarincrement(real 0.001)}
	{opt pmss}
	{opt alpha(real 0.5)}
	{opt nodraw}
	{opt nocomp:are}
	{opt col:or(string)}	
	{opt seed(int 1234)}
	{opt saving(string)}]

{pstd}
{indepvars} may contain continuous or binary variables. Where categorical variables are of interest, dummy variables for each category should be entered into {indepvars} instead.
{p_end}

{synoptset 30}{...}
{synopthdr:pmstabilityss_options}
{synoptline}
{synopt :{opt prev:alence(real)}}expected prevalence of the outcome in the development population{p_end}
{synopt :{opt cstat:istic(real 0)}}expected lower bound for the c-statistic of the developed model{p_end}
{synopt :{opt beta:s(numlist)}}specifies the expected relative weight of variables in the model{p_end}
{synopt :{opt n(numlist)}}specifies a list of additional N to investigate{p_end}
{synopt :{opt st:andardised}}standardises variables and beta effects, used instead of {opt betas()}{p_end}
{synopt :{opt lp(varname)}}linear predictor variable, used instead of {opt prev()} & {opt cstat()}{p_end}
{synopt :{opt sub:group(varname)}}specifies a categorical variable by which to investigate uncertainty{p_end}
{synopt :{opt simobs(int 500000)}}sets the N used for simulations only{p_end}
{synopt :{opt pci:width(numlist)}}specifies target interval widths within ranges of risk specified using {opt pcutpoints()}{p_end}
{synopt :{opt pcut:points(numlist)}}specifies cutpoints in the interval [0,1] corresponding with target widths specified using {opt pciwidth()}{p_end}
{synopt :{opt thresh:old(real 0)}}specifies risk threshold at which model is to be used{p_end}
{synopt :{opt tol:erance(real 0.005)}}sets tolerance for simulation achieving target prevalence and c-statistic{p_end}
{synopt :{opt logitpvarincrement(real 0.001)}}sets increment by which to iterate when identifying SE to meet the target interval widths {p_end}
{synopt :{opt pmss}}calculates uncertainty intervals at the recommended N from {cmd: pmsampsize}{p_end}
{synopt :{opt alpha(real 0.5)}}a starting value for simulation to achieve target prevalence{p_end}
{synopt :{opt nodraw}}suppress drawing of single graphs per N investigated{p_end}
{synopt :{opt nocomp:are}}suppress combining of graphs where multiple N are investigated{p_end}
{synopt :{opt col:or(string)}}set color of lines for instability plots{p_end}
{synopt :{opt seed(int 1234)}}set seed for simulations{p_end}
{synopt :{opt saving(string)}}set a suffix for saving results frames as .dta files{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}{cmd: pmstabilityss} computes the minimum sample size required for the development of a new clinical prediction model targeting precise individual risk estimates as proposed by Riley et al. 2024.
{cmd: pmstabilityss} can be used to calculate the minimum sample size for the development of models with continuous, binary or survival (time-to-event) outcomes. 
Riley et al. lay out a series of criteria the sample size should meet. These aim to minimise the overfitting and to ensure precise estimation of key parameters in the prediction model. 

{marker options}{...}



{marker examples}{...}
{title:Examples}

{phang}
We wish to develop a clinical prediction model to predict low birthweight. 
Consider that we have conducted a pilot study for our planned model development study. 
We will feed this dataset to {cmd: pmstabilityss}. 

{pstd}Using the lbw dataset as our example pilot data we can calculate some starting values for our sample size calculations. 
In particular we need an idea of the expected prevalence of our outcome of interest, a potential lower bound for how we might expect the model to perform in terms of the c-statistic, 
and some relative weights for the potential effect of each variable to be included in our new model (assuming no variable selection procedures). {p_end}
{phang2}{cmd:. }{stata webuse lbw, clear}{p_end}
{phang2}{cmd:. }{stata tab low}{p_end}
{phang2}{cmd:. }{stata logit low age lwt smoke ptl ht ui ftv}{p_end}
{phang2}{cmd:. }{stata predict pr}{p_end}
{phang2}{cmd:. }{stata roctab low pr}{p_end}

{pstd}Now we can use {cmd: pmstabilityss} to estimate the expected precision of individual risk estimates if we developed our new model on a sample of the same size as our pilot study data, 
using the following command.
Note we have set the expected relative effect of each predictor to the value seen in our pilot study.{p_end}
{phang2}{cmd:. }{stata pmstabilityss age lwt smoke ptl ht ui ftv, prev(.3122) cstat(.7382) beta(-.04328385  -.01433844    .5538123   .59443784   1.8720949   .73935229    .0235853)}{p_end}

{pstd}Given the summary table and particularly the instability plot, we would expect considerable uncertainty in individual risk estimates from a new model developed using a sample of this size (N=189). 
Let us compare this to what we would expect when using the recommended sample size using {cmd:pmsampsize}. 
To do this we can use the {opt pmss} option as below. 
{p_end}
{phang2}{cmd:. }{stata pmstabilityss age lwt smoke ptl ht ui ftv, prev(.3122) cstat(.7382) beta(-.04328385  -.01433844    .5538123   .59443784   1.8720949   .73935229    .0235853) pmss}
{p_end}

{pstd}While uncertainty in individual risk estimates is reduced somewhat by using the higher N recommended by {cmd:pmsampsize} (N=388) there is still large uncertainty, 
particularly for those with higher outcome probabilities.  
We can input an N that we believe it is feasible to obtain from our new study using the {opt n()} option. 
{p_end}
{phang2}{cmd:. }{stata pmstabilityss age lwt smoke ptl ht ui ftv, prev(.3122) cstat(.7382) beta(-.04328385  -.01433844    .5538123   .59443784   1.8720949   .73935229    .0235853) pmss n(1000)}
{p_end}

{pstd}At N=1000 there is much lower uncertainty in individual risk estimates with a median interval width of 0.08.  
We could also set a target width and allow {cmd: pmstabilityss} to calculate the required N to achieve this target.
To do this we must specify two options, first {opt pcutpoints()} and second {opt pciwidth()}. 
The following command specifies a target interval width of 0.1 for individuals with probabilities less than 0.3, 
and a target width of 0.2 for those with probabilities above 0.3. 
{p_end}
{phang2}{cmd:. }{stata pmstabilityss age lwt smoke ptl ht ui ftv, prev(.3122) cstat(.7382) beta(-.04328385  -.01433844    .5538123   .59443784   1.8720949   .73935229    .0235853) pmss pcut(0.3 1) pci(0.1 0.2)}
{p_end}

{pstd}Required sample size is significantly increased to achieve our target, but uncertainty is also dramatically reduced inline with our target.
We can also examine uncertainty in individual risk estimates across subgroups of a particular variable using the {opt subgroup()} option. 
The following command additonally asks {cmd: pmstabilityss} to investigate uncertainty by smoking status. 
{p_end}
{phang2}{cmd:. }{stata pmstabilityss age lwt smoke ptl ht ui ftv, prev(.3122) cstat(.7382) beta(-.04328385  -.01433844    .5538123   .59443784   1.8720949   .73935229    .0235853) pmss pcut(0.3 1) pci(0.1 0.2) subgroup(smoke) nodraw}
{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:pmstabilityss} stores the following in {cmd:r()} as relevant for the specific outcome type chosen:



{title:Authors}

{phang}Joie Ensor, University of Birmingham {break}
j.ensor@bham.ac.uk{p_end}

{title:Acknowledgements}

{phang}With thanks to Richard Riley{p_end}

{marker reference}{...}
{title:References}

{p 5 12 2}
Riley, R. D., Collins, G. S., Whittle, R., Archer, L., Snell, K. I., Dhiman, P., ... & Ensor, J. (2024). Sample size for developing a prediction model with a binary outcome: 
targeting precise individual risk estimates to improve clinical decisions and fairness. {it:arXiv preprint arXiv:2407.09293.}{p_end}


{title:Also see}

{psee}
Online: {helpb pmsampsize}, {helpb pmvalsampsize}, {helpb pmcalplot}
{p_end}
